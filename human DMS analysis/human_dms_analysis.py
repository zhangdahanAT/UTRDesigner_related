#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
从 DMS-Map 位点级结果中统计：
  - 每个转录本的被修饰位点数（Modified_Site_Count）
  - 每个转录本被 modified 的累计次数（Total_Modified_Reads = Σ round(pct * cov / 100)）
  - 每个转录本的长度（从 fasta 解析）
  - 被修饰的转录本总数

输出 CSV：dms_map_site_counts_with_length.csv
列：GeneFull,TranscriptID,Length,Modified_Site_Count,Total_Modified_Reads
"""

import sys
import os
import re
from collections import defaultdict

# ======== 你的默认路径（可直接运行）========
DMS_PATH_DEFAULT = r"C:\Users\Administrator\Desktop\文章各个图\supplementary figures\supplefig5\新增加的动物数据\DMS文件\GSM2241643_MegZ37_38_Mutations2.txt"
FASTA_PATH_DEFAULT = r"C:\Users\Administrator\Desktop\文章各个图\supplementary figures\supplefig5\新增加的动物数据\DMS文件\GSM2241643_human.canonical_MK.fa.txt"

# 捕获常见转录本 ID：NM_, NR_, XM_, XR_ 等
ID_PAT = re.compile(r"\b([NX][MR]_\d+\.\d+)\b")
# 备用 accession（如 V00589.1 等）
ALT_ACC_PAT = re.compile(r"\b([A-Z]{1,3}\d{5,}\.\d+)\b")

def looks_like_intish(tok: str) -> bool:
    try:
        x = float(tok)
        return abs(x - int(round(x))) < 1e-9
    except ValueError:
        return False

def looks_like_float(tok: str) -> bool:
    try:
        float(tok)
        return True
    except ValueError:
        return False

def split_fields(line: str):
    if "\t" in line:
        parts = [p for p in line.split("\t") if p != ""]
    else:
        parts = line.split()
    return parts

def pick_transcript_id(text: str) -> str:
    m = ID_PAT.search(text)
    if m:
        return m.group(1)
    m2 = ALT_ACC_PAT.search(text)
    if m2:
        return m2.group(1)
    return text

def parse_dms_map_sites(dms_path: str):
    """
    解析 DMS-Map 位点级行（第二列为 Position）。
    返回：
      site_counts: dict[GeneFull] = 位点数量
      mod_reads:   dict[GeneFull] = 累计 modified 次数（Σ round(pct * cov / 100)）
      tids:        dict[GeneFull] = TranscriptID
    兼容两种位点行格式：
      Gene  Position  PositiveMutation%  Coverage
      Gene  Position  Coverage
    """
    site_counts = defaultdict(int)
    mod_reads = defaultdict(int)
    tids = {}

    with open(dms_path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            low = line.lower()
            if low.startswith("gene"):
                # 跳过表头（包括总体统计 & 位点级表头）
                continue

            parts = split_fields(line)
            if len(parts) < 2:
                continue

            gene_full = parts[0]
            pos_token = parts[1]

            # 仅统计位点级行（第二列是 Position）
            if not looks_like_intish(pos_token):
                continue

            # 计位点数
            site_counts[gene_full] += 1
            if gene_full not in tids:
                tids[gene_full] = pick_transcript_id(gene_full)

            # 解析 Positive Mutation % 与 Coverage
            pct = None
            cov = None

            if len(parts) >= 4 and looks_like_float(parts[2]) and looks_like_float(parts[3]):
                # 典型：Gene Position Pct Coverage
                pct = float(parts[2])
                cov = float(parts[3])
            elif len(parts) >= 3 and looks_like_float(parts[2]):
                # 退化：Gene Position Coverage（没有百分比；无法算累计 modified 次数）
                cov = float(parts[2])
                pct = None
            else:
                # 无法解析，跳过累计次数计算（但位点计数保留）
                pct = None
                cov = None

            if pct is not None and cov is not None:
                # 累计 modified 次数：四舍五入为整数
                mod_reads[gene_full] += int(round(pct * cov / 100.0))

    return dict(site_counts), dict(mod_reads), tids

def parse_fasta_lengths(fasta_path: str):
    """
    解析 fasta，返回多键映射到长度的 dict：
      - header_first_token（> 后第一个空白前的 token）
      - header 中的 TranscriptID（NM_/NR_/XM_/XR_）
      - 备用 ALT_ACC_PAT（如 V00589.1）
    """
    id_to_len = {}
    header = None
    seq_len = 0

    def commit_header(h, length):
        if h is None:
            return
        first_token = h.split()[0]
        id_to_len[first_token] = length
        m = ID_PAT.search(h)
        if m:
            id_to_len[m.group(1)] = length
        m2 = ALT_ACC_PAT.search(h)
        if m2:
            id_to_len[m2.group(1)] = length

    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                # 提交上一条
                if header is not None:
                    commit_header(header, seq_len)
                header = line[1:].strip()
                seq_len = 0
            else:
                seq_len += len(line.strip())

        # 提交最后一条
        commit_header(header, seq_len)

    return id_to_len

def main():
    # 取路径：命令行优先
    if len(sys.argv) >= 3:
        dms_path = sys.argv[1]
        fasta_path = sys.argv[2]
    else:
        dms_path = DMS_PATH_DEFAULT
        fasta_path = FASTA_PATH_DEFAULT

    if not os.path.isfile(dms_path):
        print(f"[错误] 找不到 DMS 文件：{dms_path}")
        sys.exit(1)
    if not os.path.isfile(fasta_path):
        print(f"[错误] 找不到 fasta 文件：{fasta_path}")
        sys.exit(1)

    # 1) 解析 DMS 位点 & 累计 modified 次数
    site_counts, mod_reads_sum, gene_to_tid = parse_dms_map_sites(dms_path)

    # 2) 解析 fasta 长度映射
    id_to_len = parse_fasta_lengths(fasta_path)

    # 3) 合并
    rows = []
    missing_len = 0
    for gene_full, site_cnt in sorted(site_counts.items(), key=lambda x: x[1], reverse=True):
        tid = gene_to_tid.get(gene_full, gene_full)

        # 查找长度
        length = None
        if tid in id_to_len:
            length = id_to_len[tid]
        else:
            first_token = gene_full.split()[0]
            if first_token in id_to_len:
                length = id_to_len[first_token]
            else:
                alt = ALT_ACC_PAT.search(gene_full)
                if alt and alt.group(1) in id_to_len:
                    length = id_to_len[alt.group(1)]

        if length is None:
            missing_len += 1
            length = ""

        total_mod_reads = mod_reads_sum.get(gene_full, 0)

        rows.append((gene_full, tid, length, site_cnt, total_mod_reads))

    total_modified = len(site_counts)

    # 4) 打印概览
    print("=" * 90)
    print(f"被修饰的转录本总数：{total_modified}")
    print(f"其中长度未匹配到的条目数：{missing_len}")
    print("=" * 90)
    print("示例前10条：TranscriptID\tLength\tModified_Site_Count\tTotal_Modified_Reads\t(GeneFull)")
    for gene_full, tid, length, site_cnt, total_mod_reads in rows[:10]:
        print(f"{tid}\t{length}\t{site_cnt}\t{total_mod_reads}\t({gene_full})")

    # 5) 写 CSV（与 DMS 文件同目录）
    out_csv = os.path.join(os.path.dirname(os.path.abspath(dms_path)), "dms_map_site_counts_with_length.csv")
    try:
        with open(out_csv, "w", encoding="utf-8") as w:
            w.write("GeneFull,TranscriptID,Length,Modified_Site_Count,Total_Modified_Reads\n")
            for gene_full, tid, length, site_cnt, total_mod_reads in rows:
                gf = f"\"{gene_full}\"" if "," in gene_full else gene_full
                tid_csv = f"\"{tid}\"" if "," in tid else tid
                length_csv = length if isinstance(length, str) else str(length)
                w.write(f"{gf},{tid_csv},{length_csv},{site_cnt},{total_mod_reads}\n")
        print("\n已输出 CSV：", out_csv)
    except Exception as e:
        print("\n写出 CSV 失败：", e)

if __name__ == "__main__":
    main()
