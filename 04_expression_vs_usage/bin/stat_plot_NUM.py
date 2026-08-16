#!/usr/bin/env python3
"""
APA差异表达交集统计可视化脚本（每个 group1 生成独立图形）
输入: tab分隔的文本文件，包含计数数据（整数）
输出: 每个 group1 对应的 PDF 和 PNG 格式柱形图
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='为每个 group1 生成独立的 DUAPA vs DEAPA 柱形图（基于计数数据）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  python plot_apa_individual.py input.txt results/apa_plot
  python plot_apa_individual.py --width 12 --height 8 input.txt results/apa_plot
        '''
    )

    parser.add_argument('input_file', help='输入文件路径（tab分隔文本）')
    parser.add_argument('output_prefix', help='输出文件前缀（不包含扩展名，最终文件为 "前缀_group1.pdf/png"）')
    parser.add_argument('--width', type=float, default=10, help='图形宽度（英寸），默认10')
    parser.add_argument('--height', type=float, default=6, help='图形高度（英寸），默认6')
    parser.add_argument('--dpi', type=int, default=300, help='输出分辨率，默认300')
    parser.add_argument('--color_palette', default='Set2',
                       help='颜色方案，默认Set2，可选：Set2, Set3, tab20c, Pastel1等')

    return parser.parse_args()

def load_data(input_file):
    """加载数据文件"""
    try:
        df = pd.read_csv(input_file, sep='\t')
        print(f"成功加载数据，共 {len(df)} 行")
        print(f"列名: {list(df.columns)}")
        print(f"\n数据预览:")
        print(df.head())
        return df
    except Exception as e:
        print(f"加载数据失败: {e}")
        sys.exit(1)

def validate_data(df):
    """验证数据格式并确保计数列为数值类型"""
    required_cols = ['group1', 'compare', 'long∩up', 'long∩down', 'short∩up', 'short∩down']

    # 检查列名
    for col in required_cols:
        if col not in df.columns:
            print(f"错误: 缺少列 '{col}'")
            print(f"数据列名: {list(df.columns)}")
            sys.exit(1)

    # 将计数列转换为数值类型（整数）
    for col in ['long∩up', 'long∩down', 'short∩up', 'short∩down']:
        try:
            df[col] = pd.to_numeric(df[col], errors='coerce')
            if df[col].isnull().any():
                print(f"警告: 列 {col} 包含无法转换为数字的值，已替换为 NaN")
            # 检查是否所有值都是整数（可选，非必须）
            # 如果不是整数但脚本仍可运行，可以保留为浮点数
        except Exception as e:
            print(f"转换 {col} 失败: {e}")
            sys.exit(1)

    print("计数列已成功转换为数值类型")
    return df

def create_individual_plots(df, output_prefix, width=10, height=6, dpi=300, color_palette='Set2'):
    """为每个 group1 创建独立的柱形图"""

    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'Arial']
    plt.rcParams['axes.unicode_minus'] = False

    # 获取所有唯一的 group1
    groups = df['group1'].unique()
    print(f"找到 {len(groups)} 个独立组: {list(groups)}")

    # 设置颜色
    colors = plt.cm.get_cmap(color_palette, 4)

    for group in groups:
        # 获取当前组的数据
        group_df = df[df['group1'] == group].sort_values('compare')
        n_bars = len(group_df)

        if n_bars == 0:
            print(f"警告: 组 '{group}' 没有数据，跳过")
            continue

        print(f"正在绘制组: {group} (包含 {n_bars} 个比较)")

        # 创建图形
        fig, ax = plt.subplots(figsize=(width, height))

        # 设置条形图位置
        bar_width = 0.2
        index = range(n_bars)

        # 绘制四个指标
        ax.bar([i - 1.5*bar_width for i in index],
               group_df['long∩up'], bar_width,
               label='long∩up', color=colors(0), edgecolor='black')
        ax.bar([i - 0.5*bar_width for i in index],
               group_df['long∩down'], bar_width,
               label='long∩down', color=colors(1), edgecolor='black')
        ax.bar([i + 0.5*bar_width for i in index],
               group_df['short∩up'], bar_width,
               label='short∩up', color=colors(2), edgecolor='black')
        ax.bar([i + 1.5*bar_width for i in index],
               group_df['short∩down'], bar_width,
               label='short∩down', color=colors(3), edgecolor='black')

        # 设置标题和轴标签
        ax.set_title(f'DUAPA vs DEAPA - {group}', fontsize=14, fontweight='bold')
        ax.set_xlabel('compare', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)  # 改为 Count

        # 设置x轴刻度
        ax.set_xticks(index)
        ax.set_xticklabels(group_df['compare'], rotation=45, ha='right', fontsize=10)

        # 添加图例
        ax.legend(fontsize=10)

        # 添加网格
        ax.grid(True, axis='y', alpha=0.3, linestyle='--')

        # 调整y轴范围，留出空间
        y_max = group_df[['long∩up', 'long∩down', 'short∩up', 'short∩down']].max().max()
        y_max = max(y_max, 1)  # 至少1，避免y轴为0时图形为空
        ax.set_ylim(0, y_max * 1.1)

        # 调整布局
        plt.tight_layout()

        # 构建输出文件名（替换可能引起问题的字符）
        safe_group = group.replace('/', '_').replace('\\', '_').replace(' ', '_')
        pdf_file = f"{output_prefix}_{safe_group}.pdf"
        png_file = f"{output_prefix}_{safe_group}.png"

        # 保存为PDF和PNG
        plt.savefig(pdf_file, dpi=dpi, format='pdf', bbox_inches='tight')
        print(f"  保存PDF: {pdf_file}")
        plt.savefig(png_file, dpi=dpi, format='png', bbox_inches='tight')
        print(f"  保存PNG: {png_file}")

        plt.close(fig)

def main():
    """主函数"""
    args = parse_arguments()

    # 加载数据
    df = load_data(args.input_file)

    # 验证和转换数据
    df = validate_data(df)

    # 按group1和compare排序
    if 'group1' in df.columns:
        df = df.sort_values(['group1', 'compare'])

    print(f"\n数据统计:")
    print(df[['long∩up', 'long∩down', 'short∩up', 'short∩down']].describe())

    # 创建输出目录（如果需要）
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"创建输出目录: {output_dir}")

    # 生成每个组的独立图形
    print("\n开始生成各组独立图表...")
    create_individual_plots(df, args.output_prefix, args.width, args.height, args.dpi, args.color_palette)

    print(f"\n所有图表已保存，前缀为: {args.output_prefix}_[group].pdf/png")
    print("完成！")

if __name__ == "__main__":
    main()
