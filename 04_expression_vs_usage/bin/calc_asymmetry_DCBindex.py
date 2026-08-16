#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import glob
import os

# ==============================
# 参数
# ==============================
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True,
                    help="Input pattern, e.g. '*.updown_elongeshort.list'")
parser.add_argument("-o", "--output", required=True,
                    help="Output prefix")
parser.add_argument("--bootstrap", type=int, default=5000,
                    help="Bootstrap iterations (default=5000)")
args = parser.parse_args()

# ==============================
# Bootstrap for signed DCB
# ==============================
def bootstrap_signed_dcb(labels, n_boot=5000):

    if len(labels) == 0:
        return 0, 0

    boot_vals = []

    for _ in range(n_boot):
        sample = np.random.choice(labels, size=len(labels), replace=True)
        pos = np.sum(sample == 1)
        neg = np.sum(sample == -1)
        dcb = (pos - neg) / (pos + neg)
        boot_vals.append(dcb)

    lower = np.percentile(boot_vals, 2.5)
    upper = np.percentile(boot_vals, 97.5)

    return lower, upper


# ==============================
# 统计
# ==============================
results = []
files = glob.glob(args.input)

if len(files) == 0:
    print("No input files found.")
    exit()

for f in files:

    group = os.path.basename(f).split(".")[0]

    pos = 0
    neg = 0
    labels = []

    with open(f) as infile:
        for line in infile:

            line = line.strip()

            if line.startswith("组") or line.startswith("---") or line == "":
                continue

            parts = line.split()
            if len(parts) < 4:
                continue

            type_field = parts[3]

            # 正方向
            if type_field == "long-up" or type_field == "short-down":
                pos += 1
                labels.append(1)

            # 反方向
            elif type_field == "short-up" or type_field == "long-down":
                neg += 1
                labels.append(-1)

    total = pos + neg

    if total > 0:
        signed_dcb = (pos - neg) / total
        ci_low, ci_high = bootstrap_signed_dcb(np.array(labels), args.bootstrap)
    else:
        signed_dcb = 0
        ci_low, ci_high = 0, 0

    results.append([
        group,
        pos,
        neg,
        total,
        signed_dcb,
        ci_low,
        ci_high
    ])

# ==============================
# 构建表格
# ==============================
df = pd.DataFrame(results, columns=[
    "group",
    "positive_coupling",
    "negative_coupling",
    "total",
    "Signed_DCB",
    "CI_lower",
    "CI_upper"
])

# 保存完整TSV（包括Total和0组）
df.to_csv(args.output + ".tsv", sep="\t", index=False)

# ==============================
# 作图数据（去Total和total=0）
# ==============================
plot_df = df[(df["group"] != "Total") & (df["total"] > 0)].copy()

# 按 signed 值排序
plot_df = plot_df.sort_values("Signed_DCB", ascending=False)

# ==============================
# 论文级柱状图（修改后）
# ==============================
fig, ax = plt.subplots(figsize=(6, 4.5))

# 定义颜色：正值用蓝色，负值用红色（学术期刊风格）
colors = []
for val in plot_df["Signed_DCB"]:
    if val > 0:
        colors.append("#7AAAE0")   # 柔雾蓝
    else:
        colors.append("#E07B7B")   # 柔雾红

# 计算误差条的下限和上限
yerr_low = plot_df["Signed_DCB"] - plot_df["CI_lower"]
yerr_high = plot_df["CI_upper"] - plot_df["Signed_DCB"]
yerr = [yerr_low, yerr_high]

bars = ax.bar(
    plot_df["group"],
    plot_df["Signed_DCB"],
    yerr=yerr,
    capsize=4,
    color=colors,
    edgecolor="black",
    linewidth=0.6
)

# 0 中线
ax.axhline(0, color="black", linewidth=1)

# 去掉上轴线和右轴线
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# 保留左轴线和下轴线（默认保留，但为了明确可写）
ax.spines['left'].set_visible(True)
ax.spines['bottom'].set_visible(True)

ax.set_ylabel("Signed Directional Coupling Bias (DCB)")
# 去掉图标题（原代码有标题，现删除）

ax.set_xticks(range(len(plot_df)))
ax.set_xticklabels(plot_df["group"], rotation=45, ha="right")

# 自动调整 y 轴范围（基于置信区间）
ymax = max(abs(plot_df["CI_lower"]).max(),
           abs(plot_df["CI_upper"]).max())
ax.set_ylim(-ymax*1.1, ymax*1.1)

# 在每根柱子顶部打印 signed DCB 值（保留2位小数），避免与误差条重叠
offset = 0.02 * (ymax*2.2)  # 偏移量基于整个y轴范围
for i, (idx, row) in enumerate(plot_df.iterrows()):
    val = row["Signed_DCB"]
    x = i  # bar 的 x 位置索引（与 bar 顺序一致）
    # 计算上误差条端点
    upper_err = row["CI_upper"]
    # 文本放置位置：在误差条上方 offset 处，或柱子高度 + offset（取较大者）
    text_y = max(val, upper_err) + offset
    ax.text(x, text_y, f"{val:.2f}", ha='center', va='bottom', fontsize=8)

# 添加图例
positive_patch = mpatches.Patch(color="#7AAAE0", label="Elongation-associated upregulation") #Positive
negative_patch = mpatches.Patch(color="#E07B7B", label="Elongation-associated downregulation") #Negative
ax.legend(handles=[positive_patch, negative_patch], loc='upper right', frameon=False)

plt.tight_layout()

plt.savefig(args.output + ".png", dpi=600)
plt.savefig(args.output + ".pdf")
plt.close()

print("Finished.")
