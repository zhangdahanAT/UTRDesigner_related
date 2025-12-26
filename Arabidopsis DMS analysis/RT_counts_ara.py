# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 18:16:08 2025

@author: Administrator
"""

# -*- coding: utf-8 -*-
"""
统计4个样本 (2个control, 2个DMS+) 的RT次数并对齐
输出结果: rt_counts_aligned.tsv
"""

import pandas as pd
from pathlib import Path

# ========= 配置区域 =========
# 数据目录
data_dir = Path(r"C:\Users\Administrator\Desktop\2025-5-9\TAIR")

# 文件清单
files = {
    "control_551": data_dir / "sam551qian4.txt",
    "control_557": data_dir / "sam557qian4.txt",
    "dms_556":     data_dir / "sam556qian4.txt",
    "dms_552":     data_dir / "sam552qian4.txt",
}

# ========= 功能函数 =========
def load_rt_counts(fp, sample_name):
    """
    从三列文本 [FLAG, RNAME, POS] 统计 (RNAME, POS) 出现次数
    过滤条件: FLAG != 4 且 RNAME != '*' 且 POS > 0
    返回一个 Series: index = (RNAME, POS), value = count
    """
    cols = ["FLAG", "RNAME", "POS"]
    df = pd.read_csv(fp, sep=r"\s+", header=None, names=cols,
                     dtype={"FLAG": "int32", "RNAME": "string", "POS": "int64"})

    # 过滤掉未比对或无效行
    df = df[(df["FLAG"] != 4) & (df["RNAME"] != "*") & (df["POS"] > 0)]

    # 统计次数
    s = df.groupby(["RNAME", "POS"]).size()
    s.name = sample_name
    return s

# ========= 主程序 =========
def main():
    series_list = []

    for sample, fp in files.items():
        if not fp.exists():
            raise FileNotFoundError(f"找不到文件: {fp}")
        print(f"正在处理 {sample} ...")
        s = load_rt_counts(fp, sample)
        series_list.append(s)

    # 按并集对齐，缺失填0
    aligned = pd.concat(series_list, axis=1).fillna(0).astype("int64")

    # 按染色体和位置排序（字母数字混合时按字符串排序）
    aligned = aligned.sort_index(level=["RNAME", "POS"])

    # 额外加一列总数
    aligned["total"] = aligned.sum(axis=1)

    # 输出到文件
    out_path = data_dir / "rt_counts_aligned.tsv"
    aligned.to_csv(out_path, sep="\t")

    print("完成 ✅")
    print("结果保存到:", out_path)
    print("表格形状:", aligned.shape)
    print("表头:", list(aligned.columns))

if __name__ == "__main__":
    main()
    
    
    

# -*- coding: utf-8 -*-
"""
把 RT 计数与 samtools depth 对齐（只保留 1–5 号染色体）
输出: rt_depth_chr1to5.tsv
"""
import pandas as pd
from pathlib import Path

# ======= 配置 =======
data_dir = Path(r"C:\Users\Administrator\Desktop\2025-5-9\TAIR")
rt_file = data_dir / "rt_counts_aligned.tsv"

# SRR933551/557 为 control, 552/556 为 dms+（按你之前的说明）
depth_files = {
    "control_551": data_dir / "SRR933551_sorted_depth.txt",
    "control_557": data_dir / "SRR933557_sorted_depth.txt",
    "dms_556":     data_dir / "SRR933556_sorted_depth.txt",
    "dms_552":     data_dir / "SRR933552_sorted_depth.txt",
}

valid_chrs = {"1","2","3","4","5"}  # 只保留 1–5 染色体

# ======= 读取 RT（只保留 chr1-5），统一类型 =======
print("读取 RT：", rt_file)
rt_df = pd.read_csv(rt_file, sep="\t", index_col=[0, 1])
rt_df.index.set_names(["RNAME", "POS"], inplace=True)
rt_df = rt_df.reset_index()

# 统一类型 & 过滤 1-5 号染色体
rt_df["RNAME"] = rt_df["RNAME"].astype(str).str.strip()
rt_df["POS"]   = pd.to_numeric(rt_df["POS"], errors="coerce").astype("Int64")
rt_df = rt_df.dropna(subset=["POS"])
rt_df["POS"]   = rt_df["POS"].astype("int64")
rt_df = rt_df[rt_df["RNAME"].isin(valid_chrs)].copy()

print(f"RT 坐标总数（chr1-5）: {len(rt_df):,}")

# ======= 逐个 depth 文件左连接到 RT（按 RNAME, POS） =======
merged = rt_df.copy()
for sample, fp in depth_files.items():
    print(f"读取 depth：{fp}")
    d = pd.read_csv(fp, sep="\t", header=None,
                    names=["RNAME", "POS", "DEPTH"],
                    dtype={"RNAME":"string", "POS":"int64", "DEPTH":"int64"})
    # 统一类型 & 过滤 1-5 号染色体
    d["RNAME"] = d["RNAME"].astype(str).str.strip()
    d = d[d["RNAME"].isin(valid_chrs)].copy()

    # 若 depth 有重复坐标，先聚合（一般 samtools depth 不会重复，这里以防万一）
    d = d.groupby(["RNAME","POS"], as_index=False)["DEPTH"].sum()

    # 左连接（仅保留 RT 中存在的坐标）
    colname = sample + "_depth"
    merged = merged.merge(d.rename(columns={"DEPTH": colname}),
                          on=["RNAME","POS"], how="left")

    # 统计对齐命中率
    matched = merged[colname].notna().sum()
    print(f"  {sample}: 命中 {matched:,} / {len(merged):,} RT 坐标")

# 把未命中的 depth 置为 0
depth_cols = [c for c in merged.columns if c.endswith("_depth")]
merged[depth_cols] = merged[depth_cols].fillna(0).astype("int64")

# ======= 输出 =======
out_path = data_dir / "rt_depth_chr1to5.tsv"
merged = merged.set_index(["RNAME","POS"]).sort_index()
merged.to_csv(out_path, sep="\t")

print("\n完成 ✅ 输出：", out_path)
print("形状：", merged.shape)
print("列：", list(merged.columns))




## 筛选RT差值大于0的
# -*- coding: utf-8 -*-
"""
从 rt_depth_chr1to5.tsv 计算：
- DMS(556+552) 的 RT 和 depth 之和
- Control(551+557) 的 RT 和 depth 之和
- 比值差： (DMS_RT / DMS_depth) - (CTRL_RT / CTRL_depth)
仅输出差值 > 0 的行
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ======= 配置 =======
data_dir = Path(r"C:\Users\Administrator\Desktop\2025-5-9\TAIR")
in_path  = data_dir / "rt_depth_chr1to5.tsv"
out_path = data_dir / "rt_depth_chr1to5_ratio_diff_gt0.tsv"

# 期望的列名（与之前脚本一致）
rt_cols_expected = {
    "CTRL": ["control_551", "control_557"],
    "DMS":  ["dms_556", "dms_552"],
}
depth_cols_expected = {
    "CTRL": ["control_551_depth", "control_557_depth"],
    "DMS":  ["dms_556_depth", "dms_552_depth"],
}

# ======= 读取 =======
print("读取：", in_path)
df = pd.read_csv(in_path, sep="\t", index_col=[0, 1])  # MultiIndex: (RNAME, POS)

# ======= 检查列是否存在 =======
def check_columns(required_cols, df_cols, label):
    missing = [c for c in required_cols if c not in df_cols]
    if missing:
        raise KeyError(f"{label} 缺少列：{missing}\n实际列为：{list(df_cols)}")

check_columns(rt_cols_expected["CTRL"] + rt_cols_expected["DMS"], df.columns, "RT计数")
check_columns(depth_cols_expected["CTRL"] + depth_cols_expected["DMS"], df.columns, "Depth")

# ======= 求和 =======
df["dms_RT"]   = df[rt_cols_expected["DMS"]].sum(axis=1)
df["ctrl_RT"]  = df[rt_cols_expected["CTRL"]].sum(axis=1)
df["dms_depth"]  = df[depth_cols_expected["DMS"]].sum(axis=1)
df["ctrl_depth"] = df[depth_cols_expected["CTRL"]].sum(axis=1)

# ======= 计算比值（depth==0 时记为 0，避免除零）=======
def safe_ratio(numer, denom):
    # denom=0 → 0.0
    ratio = np.zeros_like(numer, dtype=float)
    mask = denom > 0
    ratio[mask] = numer[mask] / denom[mask]
    return ratio

df["ratio_dms"]  = safe_ratio(df["dms_RT"].values,  df["dms_depth"].values)
df["ratio_ctrl"] = safe_ratio(df["ctrl_RT"].values, df["ctrl_depth"].values)

# 差值 = DMS - CTRL
df["ratio_diff"] = df["ratio_dms"] - df["ratio_ctrl"]

# ======= 只保留差值 > 0 的行，并按差值降序方便查看 =======
res = df[df["ratio_diff"] > 0].copy()
res = res.sort_values("ratio_diff", ascending=False)

# ======= 只导出关键信息（如需保留原始列可去掉这一步）=======
cols_to_keep = [
    "dms_RT", "dms_depth", "ratio_dms",
    "ctrl_RT", "ctrl_depth", "ratio_ctrl",
    "ratio_diff"
]
# 若你还想保留原始4列RT和4列depth，可把上面list扩展或直接注释掉以下一行：
res = res[cols_to_keep]

# ======= 导出 =======
res.to_csv(out_path, sep="\t")
print(f"完成 ✅ 仅差值>0 的结果已保存：{out_path}")
print("行数：", len(res))
