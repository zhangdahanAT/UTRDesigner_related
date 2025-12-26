# -*- coding: utf-8 -*-
"""
Created on Wed Sep 24 09:22:44 2025

@author: Administrator
"""

# -*- coding: utf-8 -*-
"""
仅用 GTF 的 five_prime_utr 计算每个转录本的 5' 端坐标；
从 rt_depth_chr1to5.tsv 取该坐标的 dms_RT/ctrl_RT 作为 Pr(0)/Mr(0)。
——不统计无 5'UTR 的转录本——
输出: full_length_counts.tsv (transcript_id, Pr0, Mr0)
"""

import pandas as pd
from pathlib import Path

# ======== 配置 ========
data_dir = Path(r"C:\Users\Administrator\Desktop\2025-5-9\TAIR")
gtf_path = data_dir / "Arabidopsis_thaliana.TAIR10.40.gtf"
rt_path  = data_dir / "rt_depth_chr1to5.tsv"
out_path = data_dir / "full_length_counts.tsv"

valid_chrs = {"1","2","3","4","5"}
# 如果你的RT-stop坐标是“修饰位点 +1”，把下面设为 1（正链 +1、负链 −1）
STOP_OFFSET = 0

def parse_attrs(s: str):
    d={}
    for kv in s.split(";"):
        kv=kv.strip()
        if not kv:
            continue
        if " " in kv:
            k,v = kv.split(" ",1)
            d[k]=v.strip().strip('"')
    return d

# 1) 读取 GTF 中的 five_prime_utr（只保留 1-5 号染色体）
rows = []
with gtf_path.open("r", encoding="utf-8") as f:
    for line in f:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            continue
        chrom, src, feat, start, end, score, strand, frame, attrs = parts
        if feat != "five_prime_utr":
            continue
        if chrom not in valid_chrs:
            continue
        ad = parse_attrs(attrs)
        tid = ad.get("transcript_id")
        if not tid:
            continue
        rows.append((tid, chrom, int(start), int(end), strand))

utr5 = pd.DataFrame(rows, columns=["transcript_id","chrom","start1","end1","strand"])
if utr5.empty:
    raise SystemExit("GTF 中没有 five_prime_utr 条目（或被染色体过滤掉）。")

# 2) 用聚合而非 groupby.apply 取每个转录本的 5' 坐标（避免 DeprecationWarning）
#    +链取最小 start1，-链取最大 end1；同时保留 chrom 与 strand
utr5["chrom"] = utr5["chrom"].astype(str)

strand_first = utr5.groupby("transcript_id")["strand"].first()
chrom_first  = utr5.groupby("transcript_id")["chrom"].first()
pos_plus  = (utr5[utr5["strand"]== "+"].groupby("transcript_id")["start1"].min())
pos_minus = (utr5[utr5["strand"]== "-"].groupby("transcript_id")["end1"].max())
pos1 = pd.concat([pos_plus, pos_minus]).to_frame("pos1")

tx5_df = (
    pd.concat([chrom_first, strand_first, pos1], axis=1)
      .dropna(subset=["pos1"])
      .reset_index()
)

# 应用 STOP_OFFSET（正链 +offset，负链 -offset），并统一类型
tx5_df.loc[tx5_df["strand"]=="+", "pos1"] = tx5_df.loc[tx5_df["strand"]=="+", "pos1"] + STOP_OFFSET
tx5_df.loc[tx5_df["strand"]=="-", "pos1"] = tx5_df.loc[tx5_df["strand"]=="-", "pos1"] - STOP_OFFSET

tx5_df["chrom"] = tx5_df["chrom"].astype(str)
tx5_df["pos1"]  = pd.to_numeric(tx5_df["pos1"], errors="coerce").astype("Int64")
tx5_df = tx5_df.dropna(subset=["pos1"]).copy()
tx5_df["pos1"]  = tx5_df["pos1"].astype("int64")

# 3) 读取 RT-stop 表，准备 dms_RT / ctrl_RT，并统一类型（避免 merge 类型冲突）
rt_df = pd.read_csv(rt_path, sep="\t", index_col=[0,1])
rt_df.index.set_names(["chrom","pos1"], inplace=True)
rt_df = rt_df.reset_index()

rt_df["chrom"] = rt_df["chrom"].astype(str).str.strip()
rt_df["pos1"]  = pd.to_numeric(rt_df["pos1"], errors="coerce").astype("Int64")
rt_df = rt_df.dropna(subset=["pos1"]).copy()
rt_df["pos1"]  = rt_df["pos1"].astype("int64")

# 若没有现成合列就自动加和
if "dms_RT" not in rt_df.columns or "ctrl_RT" not in rt_df.columns:
    dms_cols  = [c for c in rt_df.columns if c in ("dms_556","dms_552")]
    ctrl_cols = [c for c in rt_df.columns if c in ("control_551","control_557")]
    if not dms_cols or not ctrl_cols:
        raise KeyError("rt_depth_chr1to5.tsv 未包含 dms_RT/ctrl_RT，也找不到用于求和的列")
    rt_df["dms_RT"]  = rt_df[dms_cols].sum(axis=1)
    rt_df["ctrl_RT"] = rt_df[ctrl_cols].sum(axis=1)

rt_pick = rt_df[["chrom","pos1","dms_RT","ctrl_RT"]].copy()

# 4) 左连接并把缺失置为 0（只统计有 5'UTR 的转录本）
#    这里 chrom/pos1 两边都是 str/int64，避免了你遇到的 merge 报错
fl = (
    tx5_df[["transcript_id","chrom","strand","pos1"]]
      .merge(rt_pick, on=["chrom","pos1"], how="left")
      .fillna(0)
)

fl = fl.rename(columns={"dms_RT":"Pr0","ctrl_RT":"Mr0"})[["transcript_id","Pr0","Mr0"]]
fl[["Pr0","Mr0"]] = fl[["Pr0","Mr0"]].astype("int64")

# 5) 导出
fl.to_csv(out_path, sep="\t", index=False)
print("✅ Wrote:", out_path, "rows:", len(fl))
print(fl.head(10))















# -*- coding: utf-8 -*-
"""
带进度条的 DMS reactivity 计算：
输入：
  - Arabidopsis_thaliana.TAIR10.40.gtf
  - rt_depth_chr1to5.tsv               （含 dms_RT/ctrl_RT 或复本列）
  - full_length_counts.tsv             （transcript_id, Pr0, Mr0）仅含有5'UTR转录本
输出：
  - dms_reactivity.tsv  （每个转录本的每个位点 i 的最终 DMS reactivity）
"""

import re, math
import numpy as np
import pandas as pd
from pathlib import Path


# ==== 尝试引入 tqdm（进度条）；若不可用则静默降级 ====
try:
    from tqdm.auto import tqdm
except Exception:  # 无 tqdm 则给个空实现
    def tqdm(iterable=None, **kwargs):
        return iterable if iterable is not None else []

# ===== 配置 =====
data_dir = Path(r"C:\Users\Administrator\Desktop\2025-5-9\TAIR")
gtf_path = data_dir / "Arabidopsis_thaliana.TAIR10.40.gtf"
rt_path  = data_dir / "rt_depth_chr1to5.tsv"
fl_path  = data_dir / "full_length_counts.tsv"
out_path = data_dir / "dms_reactivity.tsv"

valid_chrs = {"1","2","3","4","5"}
EPS = 1.0  # 伪计数，避免 ln(0)

# ===== 工具函数 =====
def parse_gtf_attrs(attr: str) -> dict:
    d = {}
    for m in re.finditer(r'(\S+)\s+"([^"]*)";', attr):
        d[m.group(1)] = m.group(2)
    return d

def load_exons(gtf_path: Path) -> pd.DataFrame:
    rows = []
    with gtf_path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9: continue
            chrom, src, feat, start, end, score, strand, frame, attrs = parts
            if feat != "exon": continue
            if chrom not in valid_chrs: continue
            ad = parse_gtf_attrs(attrs)
            tid = ad.get("transcript_id")
            if not tid: continue
            rows.append((tid, chrom, int(start), int(end), strand))
    ex = pd.DataFrame(rows, columns=["transcript_id","chrom","start1","end1","strand"])
    if ex.empty:
        raise SystemExit("GTF 中未解析到 exon（或被染色体过滤）。")
    ex["exon_len"] = ex["end1"] - ex["start1"] + 1
    return ex

def build_tx_model(ex: pd.DataFrame) -> pd.DataFrame:
    """
    用 GTF 的 exon 按 5'→3' 顺序拼接，得到每个转录本：
      - exon_len
      - cum_before（进入该 exon 前的累计长度）
      - tx_len（剪接后总长度 = 各 exon_len 之和）
    """
    parts = []
    for tid, g in ex.groupby("transcript_id", sort=False):
        g = g.copy()
        # 方向：正链升序，负链降序（确保是转录本 5'→3' 顺序）
        if g["strand"].iat[0] == "+":
            g = g.sort_values(["start1","end1"])
        else:
            g = g.sort_values(["start1","end1"], ascending=False)

        g["exon_len"] = g["end1"] - g["start1"] + 1
        g["tx_len"] = g["exon_len"].sum()
        # ✅ 关键修正：进入该 exon 前的累计长度
        g["cum_before"] = g["exon_len"].cumsum().shift(fill_value=0)

        parts.append(g)
    tx_exons = pd.concat(parts, ignore_index=True)
    return tx_exons


def load_rt(rt_path: Path) -> pd.DataFrame:
    df = pd.read_csv(rt_path, sep="\t", dtype={"RNAME":"string"}, index_col=[0,1])
    df.index.set_names(["chrom","pos1"], inplace=True)
    df = df.reset_index()
    df["chrom"] = df["chrom"].astype(str).str.strip()
    df["pos1"]  = pd.to_numeric(df["pos1"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["pos1"]).copy()
    df["pos1"]  = df["pos1"].astype("int64")
    df = df[df["chrom"].isin(valid_chrs)].copy()

    # 若无合列则用复本列求和
    if "dms_RT" not in df.columns or "ctrl_RT" not in df.columns:
        dms_cols  = [c for c in df.columns if c in ("dms_556","dms_552")]
        ctrl_cols = [c for c in df.columns if c in ("control_551","control_557")]
        if not dms_cols or not ctrl_cols:
            raise KeyError("rt_depth_chr1to5.tsv 缺少 dms_RT/ctrl_RT 或可加和的复本列")
        df["dms_RT"]  = df[dms_cols].sum(axis=1)
        df["ctrl_RT"] = df[ctrl_cols].sum(axis=1)

    return df[["chrom","pos1","dms_RT","ctrl_RT"]]

def map_genome_to_txpos(rt_df: pd.DataFrame, tx_ex: pd.DataFrame, keep_tids: set) -> pd.DataFrame:
    """ (chrom,pos1) 命中 exon 后计算转录本坐标 i；仅保留 keep_tids 中的转录本（带进度条） """
    tx_ex = tx_ex[tx_ex["transcript_id"].isin(keep_tids)].copy()
    if tx_ex.empty or rt_df.empty:
        return pd.DataFrame(columns=["transcript_id","i","tx_len","dms_RT","ctrl_RT"])

    # 预分组：按染色体缓存 exon 表
    ex_by_chrom = {c: g for c, g in tx_ex.groupby("chrom", sort=False)}
    out = []

    # 按染色体遍历 RT 位点，并在每个染色体内显示子进度条
    for chrom, rgrp in tqdm(rt_df.groupby("chrom", sort=False),
                            total=rt_df["chrom"].nunique(),
                            desc="坐标映射到转录本 (按染色体)"):
        eg = ex_by_chrom.get(chrom)
        if eg is None or eg.empty:
            continue
        for r in tqdm(rgrp.itertuples(index=False),
                      total=len(rgrp), leave=False, desc=f"chr{chrom}"):
            pos = int(r.pos1)
            hit = eg[(eg["start1"] <= pos) & (pos <= eg["end1"])]
            if hit.empty:
                continue
            # 命中多个 exon（多转录本/重叠）时逐一记录
            for e in hit.itertuples(index=False):
                if e.strand == "+":
                    i = int(e.cum_before + (pos - e.start1) + 1)
                else:
                    i = int(e.cum_before + (e.end1 - pos) + 1)
                out.append((e.transcript_id, i, int(e.tx_len), int(r.dms_RT), int(r.ctrl_RT)))

    if not out:
        return pd.DataFrame(columns=["transcript_id","i","tx_len","dms_RT","ctrl_RT"])

    mapped = pd.DataFrame(out, columns=["transcript_id","i","tx_len","dms_RT","ctrl_RT"])
    mapped = mapped.groupby(["transcript_id","i","tx_len"], as_index=False)[["dms_RT","ctrl_RT"]].sum()
    return mapped

def compute_reactivity(mapped: pd.DataFrame, fl_df: pd.DataFrame) -> pd.DataFrame:
    """
    Step1/2/3 计算 reactivity。只对 fl_df 中的转录本计算。
    缺失位点按 0 计入（用 ln(EPS)）；此处为向量化计算，无需进度条。
    """
    if mapped.empty:
        return mapped

    fl_df = fl_df.copy()
    fl_df["transcript_id"] = fl_df["transcript_id"].astype(str)
    mapped = mapped[mapped["transcript_id"].isin(set(fl_df["transcript_id"]))].copy()
    if mapped.empty:
        return mapped

    ln_eps = math.log(EPS)
    agg = mapped.groupby("transcript_id", as_index=False).agg(
        sum_ln_P = ("dms_RT",  lambda s: float(np.log(s.values + EPS).sum())),
        sum_ln_M = ("ctrl_RT", lambda s: float(np.log(s.values + EPS).sum())),
        n_present = ("i", "count"),
        tx_len = ("tx_len", "first")
    )
    agg = agg.merge(fl_df, on="transcript_id", how="left").fillna({"Pr0":0,"Mr0":0})

    agg["den_P"] = (np.log(agg["Pr0"] + EPS) + agg["sum_ln_P"] + (agg["tx_len"] - agg["n_present"]) * ln_eps) / agg["tx_len"]
    agg["den_M"] = (np.log(agg["Mr0"] + EPS) + agg["sum_ln_M"] + (agg["tx_len"] - agg["n_present"]) * ln_eps) / agg["tx_len"]

    denP = dict(zip(agg["transcript_id"], agg["den_P"]))
    denM = dict(zip(agg["transcript_id"], agg["den_M"]))

    mapped["P_hat"] = np.log(mapped["dms_RT"].values + EPS) / mapped["transcript_id"].map(denP).values
    mapped["M_hat"] = np.log(mapped["ctrl_RT"].values + EPS) / mapped["transcript_id"].map(denM).values
    mapped["theta"] = np.maximum(0.0, mapped["P_hat"] - mapped["M_hat"])

    # 全局 2–8% 归一化并 cap=7
    pos_theta = mapped.loc[mapped["theta"] > 0, "theta"].sort_values(ascending=False).values
    if len(pos_theta) >= 10:
        k2  = max(1, int(math.ceil(0.02 * len(pos_theta))))
        k10 = max(k2+1, int(math.ceil(0.10 * len(pos_theta))))
        mu = float(pos_theta[k2:k10].mean()) if k10 > k2 else float(pos_theta[k2-1])
    elif len(pos_theta) > 0:
        mu = float(np.median(pos_theta))
    else:
        mu = 1.0
    mapped["reactivity"] = np.minimum(7.0, mapped["theta"] / (mu if mu > 0 else 1.0))

    out = mapped[["transcript_id","i","tx_len","dms_RT","ctrl_RT","P_hat","M_hat","theta","reactivity"]].copy()
    return out.sort_values(["transcript_id","i"])

# ===== 主流程 =====
def main():
    print("读取 GTF/exon ...")
    ex = load_exons(gtf_path)

    print("构建转录本模型 ...")
    tx_ex = build_tx_model(ex)

    print("读取 full_length_counts.tsv ...")
    fl = pd.read_csv(fl_path, sep="\t")
    if not {"transcript_id","Pr0","Mr0"}.issubset(fl.columns):
        raise KeyError("full_length_counts.tsv 需要列：transcript_id, Pr0, Mr0")
    fl["transcript_id"] = fl["transcript_id"].astype(str)

    print("读取 rt_depth_chr1to5.tsv ...")
    rt = load_rt(rt_path)

    print("映射到转录本坐标 i（带进度条） ...")
    keep_tids = set(fl["transcript_id"])
    mapped = map_genome_to_txpos(rt, tx_ex, keep_tids=keep_tids)
    if mapped.empty:
        raise SystemExit("没有任何 (chr,pos) 映射到 full_length_counts 中的转录本外显子。")

    print("计算 DMS reactivity ...")
    out = compute_reactivity(mapped, fl)

    out.to_csv(out_path, sep="\t", index=False)
        
    print("✅ 完成，输出：", out_path)
    print(out.head(10))

if __name__ == "__main__":
    main()


