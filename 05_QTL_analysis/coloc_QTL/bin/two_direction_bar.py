#!/usr/bin/env python3
"""
上下相反方向柱状图（上图向上，下图向下），共享横坐标，
显示 shared/3aQTL 和 shared/3aeQTL 比例（百分比）。
下图纵轴显示整数刻度（0~max_int），柱子数字标注在底部外侧。
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import math

def main():
    parser = argparse.ArgumentParser(description='绘制上下相反方向柱状图（百分比）')
    parser.add_argument('-i', '--input', required=True, help='输入文件 (tab分隔)')
    parser.add_argument('-o', '--output', default='proportion_bar', help='输出文件前缀')
    parser.add_argument('--title', default='Proportion of shared QTLs', help='图形标题')
    parser.add_argument('--figsize', default='8,6', help='图形尺寸 (宽度,高度)')
    parser.add_argument('--dpi', type=int, default=300, help='输出分辨率')
    args = parser.parse_args()

    # 读取数据
    df = pd.read_csv(args.input, sep='\t')
    tissues = df['Tissue'].tolist()
    prop_a = (df['shared pairs'] / df["3'aQTL phenotype-Lead SNP pairs"] * 100).tolist()
    prop_ae = (df['shared pairs'] / df["3'aeQTL phenotype-Lead SNP pairs"] * 100).tolist()

    figsize = tuple(map(float, args.figsize.split(',')))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True,
                                   gridspec_kw={'height_ratios': [1, 1], 'hspace': 0})

    # 亮色
    color_a = '#FF6B6B'   # 亮红
    color_ae = '#4A9BF5'  # 亮蓝

    x = np.arange(len(tissues))
    width = 0.6

    # ---- 上图：柱子向上 ----
    bars1 = ax1.bar(x, prop_a, width, color=color_a, edgecolor='white', linewidth=0.8)
    ax1.set_ylabel('shared / 3\'aQTL (%)', fontsize=12, fontweight='bold')
    y_max_a = max(prop_a) * 1.2 if max(prop_a) > 0 else 1
    ax1.set_ylim(0, max(25, y_max_a))
    ax1.grid(False)
    ax1.set_facecolor('white')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.tick_params(bottom=False, left=True)
    ax1.yaxis.set_major_locator(plt.MultipleLocator(5))
    # 数值标签：柱顶上方
    for bar, val in zip(bars1, prop_a):
        ax1.text(bar.get_x() + bar.get_width()/2, val + 0.5,
                 f'{val:.1f}%', ha='center', va='bottom', fontsize=9)

    # ---- 下图：柱子向下，纵坐标显示正数整数刻度 ----
    neg_vals = [-v for v in prop_ae]
    bars2 = ax2.bar(x, neg_vals, width, color=color_ae, edgecolor='white', linewidth=0.8)

    max_ae = max(prop_ae)
    # 向上取整到整数，并至少为1
    max_int = max(1, math.ceil(max_ae))
    # 纵轴范围：下限为 -max_int - 0.5 以留白，上限为0
    ax2.set_ylim(-max_int - 0.5, 0)
    ax2.set_ylabel('shared / 3\'aeQTL (%)', fontsize=12, fontweight='bold')

    # 设置纵轴刻度为负整数，标签显示为正整数
    ticks_neg = [-i for i in range(0, max_int + 1)]
    ax2.set_yticks(ticks_neg)
    ax2.set_yticklabels([str(i) for i in range(0, max_int + 1)])

    ax2.grid(False)
    ax2.set_facecolor('white')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(True)
    ax2.spines['bottom'].set_visible(True)   # 保留下方脊柱作为x轴线
    ax2.tick_params(left=True, right=False, labelleft=True, labelright=False,
                    bottom=True, labelbottom=True)
    ax2.xaxis.tick_bottom()
    ax2.set_xticks(x)
    ax2.set_xticklabels(tissues, rotation=45, ha='right', fontsize=10)
    ax2.set_xlabel('Tissue', fontsize=12, fontweight='bold')

    # 数值标签：写在柱子下部顶端（即柱子底部边缘的下方）
    # bar.get_y() 是负值，其绝对值等于 val，柱子底部坐标为 -val
    # 将文字放在底部再往下偏移 0.3 个单位，确保在柱子外侧
    offset = 0.3
    for bar, val in zip(bars2, prop_ae):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_y() - offset,
                 f'{val:.1f}%', ha='center', va='top', fontsize=9)

    # 标题
    fig.suptitle(args.title, fontsize=14, fontweight='bold', y=0.95)

    plt.tight_layout()
    plt.subplots_adjust(top=0.93, bottom=0.1)
    for ext in ['.png', '.pdf']:
        outfile = f"{args.output}{ext}"
        plt.savefig(outfile, dpi=args.dpi, bbox_inches='tight')
        print(f"已保存: {outfile}")

    plt.show()

if __name__ == "__main__":
    main()
