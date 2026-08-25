import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# 1. 原始Ct数据
# ============================================================

data_raw = [
    ["ipa15′UTRcr-2", 0,   21.970273, 21.986709, 22.008405, 22.073013],
    ["ipa15′UTRcr-2", 0.5, 22.927966, 22.934163, 22.939562, 22.997605],
    ["ipa15′UTRcr-2", 1,   23.170999, 23.185058, 23.261778, 23.264669],
    ["ipa15′UTRcr-2", 2,   23.204810, 23.220850, 23.250612, 23.341058],
    ["ipa15′UTRcr-2", 4,   23.517254, 23.542194, 23.545215, 23.583161],

    ["ipa15′UTRcr-1", 0,   21.792560, 21.848504, 21.879841, 21.951487],
    ["ipa15′UTRcr-1", 0.5, 22.818439, 22.871534, 22.943239, 22.948504],
    ["ipa15′UTRcr-1", 1,   23.5354087, 23.42934986, 22.74689933, 22.84857785],
    ["ipa15′UTRcr-1", 2,   22.960363, 22.968964, 22.972368, 23.085949],
    ["ipa15′UTRcr-1", 4,   23.089128, 23.089955, 23.112555, 23.355079],

    ["NP", 0,    20.87271413, 20.95191636, 20.87441529, 20.86305979],
    ["NP", 0.5,  21.827385, 21.850916, 21.861876, 21.929965],
    ["NP", 1,    22.069421, 22.073239, 22.111508, 22.344218],
    ["NP", 2,    22.020784, 22.033153, 22.064676, 22.168898],
    ["NP", 4,    22.369215, 22.369459, 22.407631, 22.650839],

    ["ipa15′UTRin-1", 0,   21.740291, 21.791322, 21.866175, 21.881224],
    ["ipa15′UTRin-1", 0.5, 22.24198514, 22.18904316, 22.29298162, 22.21257377],
    ["ipa15′UTRin-1", 1,   23.06495007, 23.09873423, 23.02831198, 23.01869518],
    ["ipa15′UTRin-1", 2,   23.345752, 23.409182, 23.442377, 23.469525],
    ["ipa15′UTRin-1", 4,   23.432461, 23.464687, 23.516407, 23.525359],

    ["ipa15′UTRin-2", 0,   21.507843, 21.677947, 21.768227, 21.807998],
    ["ipa15′UTRin-2", 0.5, 22.451404, 22.518519, 22.609109, 22.619684],
    ["ipa15′UTRin-2", 1,   23.049956, 23.075766, 23.082743, 23.089562],
    ["ipa15′UTRin-2", 2,   22.776493, 22.787191, 22.858546, 22.897186],
    ["ipa15′UTRin-2", 4,   23.666688, 23.704310, 23.762535, 23.770694],
]

cols = ["Sample", "Time", "Ct1", "Ct2", "Ct3", "Ct4"]
df_raw = pd.DataFrame(data_raw, columns=cols)

# 固定样品和时间顺序
sample_order = ["ipa15′UTRcr-2", "ipa15′UTRcr-1", "NP", "ipa15′UTRin-1", "ipa15′UTRin-2"]
time_order = [0, 0.5, 1, 2, 4]

# ============================================================
# 2. 将宽格式Ct数据转换成长格式
# ============================================================

df_long = df_raw.melt(
    id_vars=["Sample", "Time"],
    value_vars=["Ct1", "Ct2", "Ct3", "Ct4"],
    var_name="TechnicalReplicate",
    value_name="Ct"
)

# 绝对表达量：2^(-Ct)，等价于0.5^Ct
df_long["Expression"] = 2.0 ** (-df_long["Ct"])

# ============================================================
# 3. 每个样品分别提取自己的0 min平均表达量
# ============================================================

t0_mean = (
    df_long[df_long["Time"] == 0]
    .groupby("Sample")["Expression"]
    .mean()
    .rename("T0_MeanExpression")
)

# 按样品名称匹配各自的0 min平均表达量
df_long = df_long.merge(
    t0_mean,
    left_on="Sample",
    right_index=True,
    how="left"
)

# 关键计算：
# 每个技术重复的表达量，都除以该样品自己的0 min平均表达量
df_long["RelativeExpression"] = (
    df_long["Expression"] / df_long["T0_MeanExpression"]
)

# ============================================================
# 4. 计算各样品、各时间点的平均值和标准差
# ============================================================

df_summary = (
    df_long
    .groupby(["Sample", "Time"], as_index=False)
    .agg(
        MeanCt=("Ct", "mean"),
        MeanExpression=("Expression", "mean"),
        RelativeExpression=("RelativeExpression", "mean"),
        SD=("RelativeExpression", "std"),
        N=("RelativeExpression", "count")
    )
)

# 标准误
df_summary["SEM"] = df_summary["SD"] / np.sqrt(df_summary["N"])

# 设置排序
df_summary["Sample"] = pd.Categorical(
    df_summary["Sample"],
    categories=sample_order,
    ordered=True
)

df_summary["Time"] = pd.Categorical(
    df_summary["Time"],
    categories=time_order,
    ordered=True
)

df_summary = (
    df_summary
    .sort_values(["Sample", "Time"])
    .reset_index(drop=True)
)

# ============================================================
# 5. 检查每个样品的0 min是否均归一化为1
# ============================================================

print("\n=== 每个样品自己的0 min平均表达量 ===")
print(t0_mean.to_string())

print("\n=== 0 min归一化检查（应全部等于1）===")
print(
    df_summary[df_summary["Time"] == 0][
        ["Sample", "Time", "RelativeExpression"]
    ].to_string(index=False)
)

# ============================================================
# 6. 绐制RNA decay曲线
# ============================================================

plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 1.2
plt.rcParams["xtick.major.width"] = 1.2
plt.rcParams["ytick.major.width"] = 1.2
plt.rcParams["font.size"] = 11

colors = {
    "ipa15′UTRcr-2": "#2171b5",
    "ipa15′UTRcr-1": "#fd8d3c",
    "NP":  "#31a354",
    "ipa15′UTRin-1": "#de2d26",
    "ipa15′UTRin-2": "#8856a7"
}

markers = {
    "ipa15′UTRcr-2": "o",
    "ipa15′UTRcr-1": "s",
    "NP":  "^",
    "ipa15′UTRin-1": "D",
    "ipa15′UTRin-2": "v"
}

fig, ax = plt.subplots(figsize=(7.2, 5.4), dpi=300)

for sample in sample_order:

    plot_df = df_summary[df_summary["Sample"] == sample].copy()

    # Categorical时间转换回数值，便于绘图
    x = plot_df["Time"].astype(float).to_numpy()
    y = plot_df["RelativeExpression"].to_numpy()
    yerr = plot_df["SD"].to_numpy()

    ax.errorbar(
        x,
        y,
        yerr=yerr,                 # 技术重复的标准差
        color=colors[sample],
        marker=markers[sample],
        markersize=7,
        markeredgecolor="white",
        markeredgewidth=1,
        linewidth=2,
        elinewidth=1.2,
        capsize=3,
        capthick=1.2,
        label=sample
    )

ax.set_xlabel("Time after transcription inhibition (h)", fontsize=12)
ax.set_ylabel("Relative transcript level (0 h = 1)", fontsize=12)

ax.set_xticks(time_order)
ax.set_xticklabels(["0", "0.5", "1", "2", "4"])

ax.set_xlim(-0.1, 4.1)
ax.set_ylim(bottom=0)

ax.legend(
    loc="upper right",
    frameon=False
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

# ============================================================
# 7. 保存结果
# ============================================================

plt.savefig(
    "UQB5_mRNA_decay_each_sample_T0_normalized.pdf",
    bbox_inches="tight"
)

plt.savefig(
    "UQB5_mRNA_decay_each_sample_T0_normalized.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()

# 保存数值结果
df_summary.to_csv(
    "UQB5_mRNA_decay_each_sample_T0_normalized.csv",
    index=False
)

print("\n=== 最终归一化结果 ===")
print(
    df_summary[
        [
            "Sample",
            "Time",
            "MeanCt",
            "MeanExpression",
            "RelativeExpression",
            "SD",
            "SEM"
        ]
    ].round(6).to_string(index=False)
)



#
# ============================================================
# 6. 指数衰减拟合并计算半衰期
# ============================================================

from scipy.optimize import curve_fit
from scipy.stats import t as student_t

# 一阶指数衰减模型
# 因为每个样品自己的0 min已经归一化为1，所以初始值固定为1
def decay_model(time, k):
    return np.exp(-k * time)


fit_results = []
fit_curves = {}

for sample in sample_order:

    sample_df = df_summary[
        df_summary["Sample"] == sample
    ].copy()

    x = sample_df["Time"].astype(float).to_numpy()
    y = sample_df["RelativeExpression"].astype(float).to_numpy()

    # 拟合衰减常数k
    popt, pcov = curve_fit(
        decay_model,
        x,
        y,
        p0=[0.5],
        bounds=(0, np.inf),
        maxfev=10000
    )

    k = popt[0]
    k_se = np.sqrt(pcov[0, 0])

    # 半衰期
    half_life = np.log(2) / k

    # 模型预测值
    y_pred = decay_model(x, k)

    # R²
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)

    if ss_tot > 0:
        r_squared = 1 - ss_res / ss_tot
    else:
        r_squared = np.nan

    # 95%置信区间
    # 5个时间点，模型只有1个拟合参数k
    degrees_freedom = len(x) - 1
    t_value = student_t.ppf(0.975, degrees_freedom)

    k_lower = max(k - t_value * k_se, 1e-12)
    k_upper = k + t_value * k_se

    # 由于半衰期与k呈反比，区间上下限需要反向转换
    half_life_lower = np.log(2) / k_upper
    half_life_upper = np.log(2) / k_lower

    # 用于绘制平滑拟合曲线
    x_smooth = np.linspace(
        x.min(),
        x.max(),
        300
    )

    y_smooth = decay_model(x_smooth, k)

    fit_curves[sample] = {
        "x": x_smooth,
        "y": y_smooth
    }

    fit_results.append({
        "Sample": sample,
        "Decay_constant_k": k,
        "K_SE": k_se,
        "Half_life": half_life,
        "Half_life_95CI_lower": half_life_lower,
        "Half_life_95CI_upper": half_life_upper,
        "R_squared": r_squared
    })


df_half_life = pd.DataFrame(fit_results)

print("\n=== 指数衰减拟合及半衰期结果 ===")
print(df_half_life.round(4).to_string(index=False))

# ============================================================
# 7. 绘制散点、误差线和指数拟合曲线
# ============================================================

plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 1.2
plt.rcParams["xtick.major.width"] = 1.2
plt.rcParams["ytick.major.width"] = 1.2
plt.rcParams["font.size"] = 11

colors = {
    "38d": "#2171b5",
    "45d": "#fd8d3c",
    "NP":  "#31a354",
    "in1": "#de2d26",
    "in2": "#8856a7"
}

markers = {
    "38d": "o",
    "45d": "s",
    "NP":  "^",
    "in1": "D",
    "in2": "v"
}

fig, ax = plt.subplots(figsize=(7.5, 5.7), dpi=300)

for sample in sample_order:

    plot_df = df_summary[
        df_summary["Sample"] == sample
    ].copy()

    x = plot_df["Time"].astype(float).to_numpy()
    y = plot_df["RelativeExpression"].astype(float).to_numpy()
    yerr = plot_df["SD"].astype(float).to_numpy()

    # 获取该样品的半衰期和R²
    current_result = df_half_life[
        df_half_life["Sample"] == sample
    ].iloc[0]

    half_life = current_result["Half_life"]
    r_squared = current_result["R_squared"]

    # 绘制原始时间点及技术重复标准差
    ax.errorbar(
        x,
        y,
        yerr=yerr,
        fmt=markers[sample],
        linestyle="none",
        color=colors[sample],
        markerfacecolor=colors[sample],
        markeredgecolor="white",
        markeredgewidth=1,
        markersize=7,
        elinewidth=1.2,
        capsize=3,
        capthick=1.2,
        zorder=3
    )

    # 绘制指数拟合曲线
    ax.plot(
        fit_curves[sample]["x"],
        fit_curves[sample]["y"],
        color=colors[sample],
        linewidth=2,
        label=(
            f"{sample}: "
            f"$t_{{1/2}}$ = {half_life:.2f} h, "
            f"$R^2$ = {r_squared:.2f}"
        )
    )

ax.set_xlabel(
    "Time after transcription inhibition (h)",
    fontsize=12
)

ax.set_ylabel(
    "Relative transcript level (0 h = 1)",
    fontsize=12
)

ax.set_xticks(time_order)
ax.set_xticklabels(["0", "0.5", "1", "2", "4"])

ax.set_xlim(-0.1, 4.1)
ax.set_ylim(bottom=0)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(
    loc="upper right",
    frameon=False,
    fontsize=9
)

plt.tight_layout()

# ============================================================
# 8. 保存图片和拟合结果
# ============================================================

plt.savefig(
    "UQB5_mRNA_decay_exponential_fitting.pdf",
    bbox_inches="tight"
)

plt.savefig(
    "UQB5_mRNA_decay_exponential_fitting.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()

# 保存半衰期结果
df_half_life.to_csv(
    "UQB5_mRNA_half_life_results.csv",
    index=False
)

# 保存归一化表达量
df_summary.to_csv(
    "UQB5_mRNA_decay_normalized_data.csv",
    index=False
)

# Sample  Decay_constant_k   K_SE  Half_life  Half_life_95CI_lower  Half_life_95CI_upper  R_squared
#    38d            0.5666 0.1827     1.2233                0.6455          1.167660e+01     0.4912
#    45d            0.5106 0.1955     1.3576                0.6581          6.931472e+11     0.2056
#     NP            0.5693 0.1914     1.2175                0.6297          1.828530e+01     0.4494
#    in1            0.5526 0.1211     1.2544                0.7799          3.204200e+00     0.8199
#    in2            0.6255 0.1794     1.1081                0.6169          5.439400e+00     0.6678

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit, least_squares
from scipy.stats import f

# ============================================================
# 1. 将宽格式Ct数据转成长格式
# ============================================================

df_long = df_raw.melt(
    id_vars=["Sample", "Time"],
    value_vars=["Ct1", "Ct2", "Ct3", "Ct4"],
    var_name="Replicate",
    value_name="Ct"
)

# ============================================================
# 2. 每个样品、每个重复分别以自己的0 h为基准
#    Relative abundance = 2^-(Ct_t - Ct_0)
# ============================================================

df_t0 = (
    df_long[df_long["Time"] == 0]
    [["Sample", "Replicate", "Ct"]]
    .rename(columns={"Ct": "Ct0"})
)

df_long = df_long.merge(
    df_t0,
    on=["Sample", "Replicate"],
    how="left"
)

df_long["DeltaCt"] = df_long["Ct"] - df_long["Ct0"]

df_long["Relative_abundance"] = 2 ** (-df_long["DeltaCt"])


# ============================================================
# 3. 定义one-phase exponential decay模型
#
#    Y = exp(-k*t)
#
#    因为0 h已经归一化为1，因此Y0固定为1
# ============================================================

def decay_model(t, k):
    return np.exp(-k * t)


# ============================================================
# 4. 每个材料单独拟合k和half-life
# ============================================================

fit_results = []

for sample in sample_order:

    temp = df_long[df_long["Sample"] == sample]

    x = temp["Time"].values.astype(float)
    y = temp["Relative_abundance"].values.astype(float)

    popt, pcov = curve_fit(
        decay_model,
        x,
        y,
        p0=[0.5],
        bounds=(0, np.inf)
    )

    k = popt[0]

    # k标准误
    k_se = np.sqrt(np.diag(pcov))[0]

    # half-life
    half_life = np.log(2) / k

    # 95% CI of k
    k_lower = k - 1.96 * k_se
    k_upper = k + 1.96 * k_se

    fit_results.append({
        "Sample": sample,
        "k": k,
        "k_SE": k_se,
        "k_95CI_lower": k_lower,
        "k_95CI_upper": k_upper,
        "Half_life_h": half_life
    })


fit_df = pd.DataFrame(fit_results)

print("\n========================================")
print("Individual decay fitting")
print("========================================")

print(
    fit_df.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)


# ============================================================
# 5. 辅助函数：单个材料拟合k并计算RSS
# ============================================================

def fit_single_group(group):

    t = group["Time"].values.astype(float)
    y = group["Relative_abundance"].values.astype(float)

    result = least_squares(
        lambda p: y - np.exp(-p[0] * t),
        x0=[0.5],
        bounds=(0, np.inf)
    )

    k = result.x[0]

    RSS = np.sum(result.fun ** 2)

    return k, RSS


# ============================================================
# 6. Extra sum-of-squares F test
#
#    Restricted model:
#        WT和mutant共用一个k
#
#    Full model:
#        WT和mutant分别拥有自己的k
#
#    H0:
#        k_WT = k_mutant
#
# ============================================================

def compare_k_extra_sum_squares(
    dataframe,
    control,
    mutant
):

    combined = dataframe[
        dataframe["Sample"].isin([control, mutant])
    ].copy()

    # --------------------------------------------------------
    # Restricted model:
    # 两个材料共用同一个k
    # --------------------------------------------------------

    t_all = combined["Time"].values.astype(float)

    y_all = combined[
        "Relative_abundance"
    ].values.astype(float)

    restricted_fit = least_squares(
        lambda p: y_all - np.exp(-p[0] * t_all),
        x0=[0.5],
        bounds=(0, np.inf)
    )

    k_shared = restricted_fit.x[0]

    RSS_restricted = np.sum(
        restricted_fit.fun ** 2
    )


    # --------------------------------------------------------
    # Full model:
    # control和mutant分别拟合自己的k
    # --------------------------------------------------------

    control_data = combined[
        combined["Sample"] == control
    ]

    mutant_data = combined[
        combined["Sample"] == mutant
    ]

    k_control, RSS_control = fit_single_group(
        control_data
    )

    k_mutant, RSS_mutant = fit_single_group(
        mutant_data
    )

    RSS_full = RSS_control + RSS_mutant


    # --------------------------------------------------------
    # Extra sum-of-squares F test
    # --------------------------------------------------------

    N = len(combined)

    # restricted model只有1个参数k
    p_restricted = 1

    # full model有2个参数：
    # k_control + k_mutant
    p_full = 2

    df_num = p_full - p_restricted

    df_den = N - p_full

    F_value = (
        (RSS_restricted - RSS_full)
        / df_num
    ) / (
        RSS_full / df_den
    )

    P_value = f.sf(
        F_value,
        df_num,
        df_den
    )

    return {
        "Comparison": f"{mutant} vs {control}",
        "k_control": k_control,
        "k_mutant": k_mutant,
        "k_shared": k_shared,
        "F": F_value,
        "df1": df_num,
        "df2": df_den,
        "P_value": P_value
    }


# ============================================================
# 7. 四个突变体分别与NP比较
# ============================================================

control = "NP"

mutants = [
    "ipa15′UTRcr-2",
    "ipa15′UTRcr-1",
    "ipa15′UTRin-1",
    "ipa15′UTRin-2"
]

comparison_results = []

for mutant in mutants:

    result = compare_k_extra_sum_squares(
        df_long,
        control,
        mutant
    )

    comparison_results.append(result)


comparison_df = pd.DataFrame(
    comparison_results
)


# ============================================================
# 8. 添加显著性判断
# ============================================================

comparison_df["Significance"] = np.where(
    comparison_df["P_value"] < 0.05,
    "Significant",
    "ns"
)


print("\n========================================")
print("Extra sum-of-squares F test")
print("========================================")

print(
    comparison_df.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)


# ============================================================
# 9. 输出Excel
# ============================================================

with pd.ExcelWriter(
    "IPA1_mRNA_decay_statistics.xlsx"
) as writer:

    fit_df.to_excel(
        writer,
        sheet_name="Decay_parameters",
        index=False
    )

    comparison_df.to_excel(
        writer,
        sheet_name="F_test",
        index=False
    )

    df_long.to_excel(
        writer,
        sheet_name="Normalized_data",
        index=False
    )


print(
    "\nOutput: IPA1_mRNA_decay_statistics.xlsx"
)









