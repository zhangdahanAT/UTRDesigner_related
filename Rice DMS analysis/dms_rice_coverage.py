# -*- coding: utf-8 -*-
"""
Created on Sun Aug 17 16:47:00 2025

@author: Administrator
"""

# -*- coding: utf-8 -*-
"""
Rice 5′UTR read count 统计脚本（IntervalTree加速，基因交并集）
- 读取水稻GFF/GTF（自动识别 *.gff*）
- 统计四个FLAG/RNAME/POS文本中read起点落在5′UTR的次数（命中多个片段全部计入）
- 归一化为 hits/UTR_length
- 以阈值 >1 统计每文件阳性转录本数，并对 SRR6010248 & SRR6010249 计算【基因】的交/并集大小
"""

import os, re, glob, time
import numpy as np
import pandas as pd
from collections import defaultdict
from intervaltree import IntervalTree
from tqdm import tqdm

# ===================== 配置 =====================
root_dir = r"C:\Users\Administrator\Desktop\2025-5-9\rice"

sam_files = [
    os.path.join(root_dir, "SRR6010248_2to4lie.txt"),
    os.path.join(root_dir, "SRR6010249_2to4lie.txt"),
    os.path.join(root_dir, "SRR6010250_2to4lie.txt"),
    os.path.join(root_dir, "SRR6010251_2to4lie.txt"),
]
# 自动选取目录下第一个 *.gff* 作为注释（如需固定，可写死路径）
gff_candidates = glob.glob(os.path.join(root_dir, "*.gff*"))
if not gff_candidates:
    raise FileNotFoundError(f"在目录中未找到 GFF/GTF：{root_dir}\\*.gff*")
gff_file = gff_candidates[0]
print(f"[INFO] 使用注释文件：{gff_file}")

output_excel = os.path.join(root_dir, "5UTR_rice_counts.xlsx")
THRESHOLD = 1.0  # 用于 hits/UTR_length 的阈值
PAIR_FOR_SET = ("SRR6010248_2to4lie.txt", "SRR6010249_2to4lie.txt")  # 用于基因交并集的两个文件

# ===================== 工具函数 =====================
def norm_chrom(name: str) -> str:
    """归一化染色体名：chr01/Chr01/01 -> 1；保留 Mt/Pt 等"""
    if not name:
        return name
    s = name.strip()
    s = re.sub(r'^(chr|Chr)', '', s)
    # 若是纯数字串（可能带前导0），去掉前导0
    if re.fullmatch(r'\d+', s):
        s = re.sub(r'^0+', '', s)
        if s == '':
            s = '0'
    return s

def parse_transcript_from_attr(attr: str):
    """从 five_prime_UTR 的属性解析转录本ID：支持 Parent=xxx 或 transcript_id(=| " ")"""
    if not attr:
        return None
    m = re.search(r'\btranscript_id\s*=\s*([^;]+)', attr)
    if m:
        return m.group(1).strip()
    m = re.search(r'\btranscript_id\s+"([^"]+)"', attr)
    if m:
        return m.group(1).strip()
    m = re.search(r'\bParent=([^;]+)', attr)
    if m:
        parent = m.group(1).split(',')[0].strip()
        parent = re.sub(r'^(transcript:|mRNA:)', '', parent)
        return parent
    return None

def parse_mrna_tx_gene(attr: str):
    """从 mRNA 的属性里解析 转录本ID 与 基因ID（优先 Locus_id=Os01g...；次选 gene_id/gene=）"""
    if not attr:
        return None, None
    # 转录本ID（mRNA的ID）
    tx = None
    m = re.search(r'\bID=([^;]+)', attr)
    if m:
        tx = m.group(1).strip()
    # 基因ID
    gene = None
    m = re.search(r'\bLocus_id=([^;]+)', attr)  # 例如 Os01g0100100
    if m:
        gene = m.group(1).strip()
    else:
        m = re.search(r'\bgene_id=([^;]+)', attr)
        if m:
            gene = m.group(1).strip()
        else:
            m = re.search(r'\bgene=([^;]+)', attr)
            if m:
                gene = m.group(1).strip()
    return tx, gene

def count_lines(path: str) -> int:
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        return sum(1 for _ in f)

# ===================== 1) 解析注释：5′UTR + 转录本→基因映射 =====================
utr_records = []  # (chrom_norm, start, end, strand, transcript_id)
tx2gene = {}      # {transcript_id: gene_id}

with open(gff_file, "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        chrom, source, ftype, start, end, score, strand, phase, attr = parts

        # 收集 mRNA 的 transcript→gene
        if ftype.lower() == "mrna":
            tx, gene = parse_mrna_tx_gene(attr)
            if tx:
                tx2gene.setdefault(tx, gene)
            continue

        # 收集 five_prime_UTR 区间
        if ftype.lower() == "five_prime_utr":
            tid = parse_transcript_from_attr(attr)
            if not tid:
                continue
            try:
                start_i = int(start); end_i = int(end)
            except ValueError:
                continue
            utr_records.append((norm_chrom(chrom), start_i, end_i, strand, tid))

if not utr_records:
    raise RuntimeError("GFF/GTF 中未解析到 five_prime_UTR 记录，请检查注释文件。")

utr_df = pd.DataFrame(utr_records, columns=["chrom", "start", "end", "strand", "transcript_id"])
utr_df["length"] = utr_df["end"] - utr_df["start"] + 1

# 每转录本 UTR 总长度
utr_len = (
    utr_df.groupby("transcript_id", as_index=False)["length"]
    .sum()
    .rename(columns={"length": "5UTR_length"})
)
# 加上 gene_id 列（若无可为空）
utr_len["gene_id"] = utr_len["transcript_id"].map(tx2gene)

# ===================== 2) 构建 IntervalTree（按染色体） =====================
trees = defaultdict(IntervalTree)
for _, row in utr_df.iterrows():
    trees[row["chrom"]].addi(int(row["start"]), int(row["end"]) + 1, row["transcript_id"])

print(f"[INFO] 构建IntervalTree：{len(trees)}条染色体，{len(utr_df)}个UTR片段，{len(utr_len)}个转录本")

# ===================== 3) 逐文件统计 =====================
result_df = utr_len.copy()

for idx, sam_file in enumerate(sam_files, start=1):
    t0 = time.time()
    counts = defaultdict(int)

    total = count_lines(sam_file)
    with open(sam_file, "r", encoding="utf-8", errors="ignore") as f, \
         tqdm(total=total, desc=f"Scanning {os.path.basename(sam_file)}", unit="line") as pbar:

        for line in f:
            pbar.update(1)
            s = line.strip()
            if not s:
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            flag, rname, pos = parts[0], parts[1], parts[2]

            # 跳过未比对/无位点
            if flag == "4" or rname == "*":
                continue
            try:
                pos_i = int(pos)
            except ValueError:
                continue

            tree = trees.get(norm_chrom(rname))
            if not tree:
                continue

            # 命中多个片段全部计入
            for iv in tree.at(pos_i):
                counts[iv.data] += 1

    col = f"count_file_{idx}"
    result_df[col] = result_df["transcript_id"].map(lambda x: counts.get(x, 0)).astype(int)

    dt = time.time() - t0
    print(f"[DONE] {os.path.basename(sam_file)} in {dt:.2f}s; "
          f"nonzero transcripts: {(result_df[col] > 0).sum()}")

# ===================== 4) 归一化 & 阈值统计（> THRESHOLD） =====================
len_safe = result_df["5UTR_length"].replace(0, np.nan)

norm_cols = []
gt_tx_sets = {}   # {basename: set(transcript_id)}
gt_gene_sets = {} # {basename: set(gene_id)}

for i, sam in enumerate(sam_files, start=1):
    c = f"count_file_{i}"
    nc = f"norm_file_{i}"
    result_df[nc] = result_df[c] / len_safe
    norm_cols.append(nc)

    base = os.path.basename(sam)
    mask = result_df[nc] > THRESHOLD

    # 满足阈值的转录本集合
    tids = set(result_df.loc[mask, "transcript_id"])
    gt_tx_sets[base] = tids

    # 满足阈值的基因集合（去掉缺失gene_id）
    gids = set(result_df.loc[mask, "gene_id"].dropna().astype(str))
    gt_gene_sets[base] = gids

    print(f"{base}: hits/UTR_length > {THRESHOLD} 的转录本数 = {len(tids)}；基因数 = {len(gids)}")

# 仅对 SRR6010248 与 SRR6010249 计算交/并集（基因）
A_name, B_name = PAIR_FOR_SET
if A_name not in gt_gene_sets or B_name not in gt_gene_sets:
    raise KeyError(f"PAIR_FOR_SET 中的文件名未在 sam_files 中找到：{PAIR_FOR_SET}")

gene_inter = gt_gene_sets[A_name] & gt_gene_sets[B_name]
gene_union = gt_gene_sets[A_name] | gt_gene_sets[B_name]

# 同时给出转录本交并集，便于核对
tx_inter = gt_tx_sets[A_name] & gt_tx_sets[B_name]
tx_union = gt_tx_sets[A_name] | gt_tx_sets[B_name]

print(f"\n[汇总-基因] {A_name} 与 {B_name} 交集(>{THRESHOLD}) 基因数 = {len(gene_inter)}")
print(f"[汇总-基因] {A_name} 与 {B_name} 并集(>{THRESHOLD}) 基因数 = {len(gene_union)}")
print(f"[汇总-转录本] 交集 = {len(tx_inter)}；并集 = {len(tx_union)}")

# ===================== 5) 输出 Excel =====================
# 主表：Transcript_ID, gene_id, 5UTR_length, counts + norms
cols_main = ["transcript_id", "gene_id", "5UTR_length",
             "count_file_1","count_file_2","count_file_3","count_file_4"] + norm_cols
out_df = result_df[cols_main].copy()
out_df.rename(columns={"transcript_id":"Transcript_ID"}, inplace=True)

with pd.ExcelWriter(output_excel, engine="openpyxl") as w:
    out_df.to_excel(w, sheet_name="counts", index=False)

    # 交/并集清单——基因
    pd.DataFrame(sorted(gene_inter), columns=["gene_id"]).to_excel(
        w, sheet_name=f"{A_name}&{B_name}_genes_inter_gt{THRESHOLD}", index=False
    )
    pd.DataFrame(sorted(gene_union), columns=["gene_id"]).to_excel(
        w, sheet_name=f"{A_name}|{B_name}_genes_union_gt{THRESHOLD}", index=False
    )

    # 交/并集清单——转录本（可选）
    pd.DataFrame(sorted(tx_inter), columns=["Transcript_ID"]).to_excel(
        w, sheet_name=f"{A_name}&{B_name}_tx_inter_gt{THRESHOLD}", index=False
    )
    pd.DataFrame(sorted(tx_union), columns=["Transcript_ID"]).to_excel(
        w, sheet_name=f"{A_name}|{B_name}_tx_union_gt{THRESHOLD}", index=False
    )

print(f"\n结果已保存：{output_excel}")
