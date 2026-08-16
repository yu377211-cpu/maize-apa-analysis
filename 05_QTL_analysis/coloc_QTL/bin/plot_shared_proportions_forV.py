#!/usr/bin/env python3
"""
根据输入文件（制表符分隔）绘制两张图：
  列顺序：Tissue, A_total, B_total, AB_overlap
  1) 比例堆叠柱状图（每个组织一根柱子，总高100%，显示 A-only, shared, B-only 比例）
     union = A + B - AB
  2) Shared比例分组柱状图（比较 shared/A 和 shared/B）

输出 PNG 和 PDF。
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(description='绘制 shared 比例图（通用列顺序）')
    parser.add_argument('-i', '--input', required=True, help='输入文件 (制表符分隔，列顺序: Tissue, A_total, B_total, AB_overlap)')
    parser.add_argument('-o', '--output', default='shared_proportion', help='输出文件前缀')
    parser.add_argument('--figsize', default='10,6', help='图形尺寸 (宽度,高度)')
    parser.add_argument('--dpi', type=int, default=300, help='输出分辨率')
    args = parser.parse_args()

    # 读取文件，不依赖列名，直接按列位置读取
    df = pd.read_csv(args.input, sep='\t', header=0)
    # 假设列顺序：Tissue, A_total, B_total, AB_overlap
    tissues = df.iloc[:, 0].tolist()
    A = df.iloc[:, 1]
    B = df.iloc[:, 2]
    AB = df.iloc[:, 3]

    # 计算各组成部分
    a_only = A - AB
    b_only = B - AB
    shared = AB
    union = A + B - AB   # 总唯一数

    # ---- 图1：比例堆叠柱状图 ----
    fig1, ax1 = plt.subplots(figsize=tuple(map(float, args.figsize.split(','))))

    a_pct = a_only / union * 100
    shared_pct = shared / union * 100
    b_pct = b_only / union * 100

    #color_a_only = '#88C179'   # 绿色
    #color_shared = '#02B3EA'   # 蓝色
    #color_b_only = '#FBE729'   # 黄色
    color_a_only = '#C2FFAD'
    color_shared = '#94EAFF'
    color_b_only = '#FFAFAF'

    x = np.arange(len(tissues))
    width = 0.6

    # 堆叠：a-only -> shared -> b-only
    ax1.bar(x, a_pct, width, label="3aQTL-only", color=color_a_only)
    ax1.bar(x, shared_pct, width, bottom=a_pct, label='Shared', color=color_shared)
    ax1.bar(x, b_pct, width, bottom=a_pct + shared_pct, label="3aeQTL-only", color=color_b_only)

    # 添加百分比标签
    for i in range(len(tissues)):
        ax1.text(x[i], a_pct[i]/2, f'{a_pct[i]:.1f}%', ha='center', va='center', fontsize=9, color='black')
        y_shared_mid = a_pct[i] + shared_pct[i]/2
        ax1.text(x[i], y_shared_mid, f'{shared_pct[i]:.1f}%', ha='center', va='center', fontsize=9, color='black')
        y_b_mid = a_pct[i] + shared_pct[i] + b_pct[i]/2
        ax1.text(x[i], y_b_mid, f'{b_pct[i]:.1f}%', ha='center', va='center', fontsize=9, color='black')

    ax1.set_xticks(x)
    ax1.set_xticklabels(tissues)
    ax1.set_xlabel('Tissue', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Percentage of total unique variants (%)', fontsize=12, fontweight='bold')
    ax1.set_ylim(0, 100)
    ax1.set_title('Composition of shared and unique variants', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper left', framealpha=0.9)

    ax1.set_facecolor('white')
    ax1.grid(False)

    plt.tight_layout()
    fig1.savefig(f"{args.output}_stacked.png", dpi=args.dpi, bbox_inches='tight')
    fig1.savefig(f"{args.output}_stacked.pdf", bbox_inches='tight')
    print(f"保存堆叠图: {args.output}_stacked.png/pdf")

    # ---- 图2：Shared比例分组柱状图 ----
    fig2, ax2 = plt.subplots(figsize=tuple(map(float, args.figsize.split(','))))

    prop_a = shared / A * 100
    prop_b = shared / B * 100

    x2 = np.arange(len(tissues))
    width2 = 0.35

    bars_a = ax2.bar(x2 - width2/2, prop_a, width2, label='shared / 3aQTL', color=color_a_only)
    bars_b = ax2.bar(x2 + width2/2, prop_b, width2, label='shared / 3aeQTL', color=color_b_only)

    for bar, val in zip(bars_a, prop_a):
        ax2.text(bar.get_x() + bar.get_width()/2, val + 0.5, f'{val:.1f}%', ha='center', va='bottom', fontsize=8)
    for bar, val in zip(bars_b, prop_b):
        ax2.text(bar.get_x() + bar.get_width()/2, val + 0.5, f'{val:.1f}%', ha='center', va='bottom', fontsize=8)

    ax2.set_xticks(x2)
    ax2.set_xticklabels(tissues)
    ax2.set_xlabel('Tissue', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Shared variants proportion (%)', fontsize=12, fontweight='bold')
    ax2.set_title('Comparison of shared proportion', fontsize=14, fontweight='bold')
    ax2.legend(loc='upper left', framealpha=0.9)

    ax2.set_facecolor('white')
    ax2.grid(False)

    plt.tight_layout()
    fig2.savefig(f"{args.output}_proportion.png", dpi=args.dpi, bbox_inches='tight')
    fig2.savefig(f"{args.output}_proportion.pdf", bbox_inches='tight')
    print(f"保存比例图: {args.output}_proportion.png/pdf")

    print("完成！")

if __name__ == "__main__":
    main()
