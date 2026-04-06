# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 09:30:12 2025

@author: Administrator
"""

# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, pearsonr
from matplotlib.colors import LinearSegmentedColormap, Normalize

# ---------- 工具 ----------
def normalize_cols(cols):
    out = []
    for c in cols:
        if c is None:
            out.append(c)
            continue
        s = str(c).strip().lower()
        s = s.replace('\ufeff', '')            # 去BOM
        s = s.replace(' ', '_')                # 空格->下划线
        s = s.replace('__', '_')
        out.append(s)
    return out

def find_col(df, targets):
    """在 df.columns 里按候选名列表匹配第一命中的列名；支持 startswith/完全相等"""
    cols = list(df.columns)
    # 先完全相等
    for t in targets:
        for c in cols:
            if c == t:
                return c
    # 再 startswith
    for t in targets:
        for c in cols:
            if c.startswith(t):
                return c
    raise KeyError(f"找不到目标列（候选：{targets}）; 当前列名：{cols[:10]} ...")

def gc_at_stems(seq: str, struct: str) -> float:
    if not isinstance(seq, str) or not isinstance(struct, str):
        return np.nan
    n = min(len(seq), len(struct))
    if n == 0:
        return np.nan
    seq_u = seq[:n].upper()
    struct_s = struct[:n]
    mask = [(ch == '(') or (ch == ')') for ch in struct_s]
    if not any(mask):
        return np.nan
    bases = [b for b, m in zip(seq_u, mask) if m]
    if not bases:
        return np.nan
    gc = sum(b in ('G', 'C') for b in bases)
    return gc / len(bases)

def density_scatter(x, y, title, xlab, ylab, pdf_out):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # 仅保留有限且位于 [0,1] 的点（含边界）
    m = (
        np.isfinite(x) & np.isfinite(y) &
        (x >= 0.0) & (x <= 1.0) &
        (y >= 0.0) & (y <= 1.0)
    )
    x, y = x[m], y[m]
    if len(x) == 0:
        print("No valid points in [0,1] to plot.")
        return

    # KDE 密度
    xy = np.vstack([x, y])
    kde = gaussian_kde(xy)
    dens = kde(xy)

    cmap = LinearSegmentedColormap.from_list(
        "gray_blue_red",
        [(0.8, 0.8, 0.8, 0.3), (0.1, 0.1, 1), (1, 0, 0)],
        N=256
    )
    vmin = np.quantile(dens, 0.01)
    vmax = np.quantile(dens, 0.99)
    norm = Normalize(vmin=vmin, vmax=vmax)

    plt.figure(figsize=(8, 8))
    sc = plt.scatter(x, y, c=dens, cmap=cmap, norm=norm, s=20, edgecolors='none')
    plt.colorbar(sc, label="Density")
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.axis('square')

    # 标注 PCC & n
    pcc, pval = pearsonr(x, y)
    plt.text(0.05, 0.95, f"PCC = {pcc:.3f}\nn = {len(x)}",
             transform=plt.gca().transAxes, ha='left', va='top', fontsize=12)

    plt.tight_layout()
    plt.savefig(pdf_out, format='pdf', bbox_inches='tight')


# ---------- 读取 Excel，仅取所需四列 ----------
rice_xlsx = r"C:\Users\Administrator\Desktop\2025-5-9\rice\aligned_with_INVIVOstructures.xlsx"
df_raw = pd.read_excel(rice_xlsx, sheet_name=0, dtype=str)   # 以字符串读，避免丢字符
df_raw.columns = normalize_cols(df_raw.columns)

# 可能的列名别名（做容错）
seq_col   = find_col(df_raw, ["sequence"])
stru_col  = find_col(df_raw, ["structure"])
seqi_col  = find_col(df_raw, ["sequence_invivo", "sequence_in_vivo", "sequence_in-vivo"])
strui_col = find_col(df_raw, ["structure_invivo", "structure_in_vivo", "structure_in-vivo"])

invivo_merge = df_raw[[seq_col, stru_col, seqi_col, strui_col]].rename(columns={
    seq_col:   "sequence",
    stru_col:  "structure",
    seqi_col:  "sequence_invivo",
    strui_col: "structure_invivo",
}).copy()

# 统一大写序列
invivo_merge["sequence"]        = invivo_merge["sequence"].str.upper()
invivo_merge["sequence_invivo"] = invivo_merge["sequence_invivo"].str.upper()

# ---------- 计算茎区 GC 含量 ----------
invivo_merge["gc_stem_UTR5"]   = invivo_merge.apply(lambda r: gc_at_stems(r["sequence"],        r["structure"]),         axis=1)
invivo_merge["gc_stem_INVIVO"] = invivo_merge.apply(lambda r: gc_at_stems(r["sequence_invivo"], r["structure_invivo"]), axis=1)

# 清洗并导出
invivo_merge = invivo_merge.replace([np.inf, -np.inf], np.nan).dropna(subset=["gc_stem_UTR5", "gc_stem_INVIVO"])
invivo_merge.to_csv("invivo_merge_rice_stem_GC_fraction.csv", index=False)

# ---------- 打印 PCC ----------
pcc, pval = pearsonr(invivo_merge["gc_stem_UTR5"], invivo_merge["gc_stem_INVIVO"])
print(f"PCC = {pcc:.3f}, P-value = {pval:.3e}, n = {len(invivo_merge)}")

# ---------- 作图 ----------
density_scatter(
    invivo_merge["gc_stem_UTR5"],
    invivo_merge["gc_stem_INVIVO"],
    "GC fraction at stem positions: UTR5 vs In vivo (rice)",
    "UTR5 stem GC fraction",
    "In vivo stem GC fraction",
    "UTR5_vs_INVIVO_stem_GC_scatter_rice.pdf"
)

plt.show()
