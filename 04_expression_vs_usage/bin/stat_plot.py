#!/usr/bin/env python3
"""
APA差异表达交集统计可视化脚本
输入: tab分隔的文本文件，包含百分比数据
输出: PDF和PNG格式的柱形图
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='sig_DUAPA vs sig_DEAPA',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  python plot_apa_stats.py input.txt output_prefix
  python plot_apa_stats.py --width 12 --height 8 input.txt results/apa_plot
        '''
    )
    
    parser.add_argument('input_file', help='输入文件路径（tab分隔文本）')
    parser.add_argument('output_prefix', help='输出文件前缀（不包含扩展名）')
    parser.add_argument('--width', type=float, default=14, help='图形宽度（英寸），默认14')
    parser.add_argument('--height', type=float, default=8, help='图形高度（英寸），默认8')
    parser.add_argument('--dpi', type=int, default=300, help='输出分辨率，默认300')
    parser.add_argument('--color_palette', default='Set2', 
                       help='颜色方案，默认Set2，可选：Set2, Set3, tab20c, Pastel1等')
    parser.add_argument('--sort_by', default='group1', 
                       help='排序依据，默认group1，可选：group1, compare')
    parser.add_argument('--style', default='single', 
                       help='图表样式，默认single，可选：single, grouped, stacked, all')
    
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
    """验证数据格式"""
    required_cols = ['group1', 'compare', 'long∩up', 'long∩down', 'short∩up', 'short∩down']
   #long∩up	long∩down	short∩up	short∩down 
    # 检查列名
    for col in required_cols:
        if col not in df.columns:
            print(f"错误: 缺少列 '{col}'")
            print(f"数据列名: {list(df.columns)}")
            sys.exit(1)
    
    # 转换百分比为数值
    for col in ['long∩up', 'long∩down', 'short∩up', 'short∩down']:
        try:
            # 移除百分号并转换为浮点数
            df[col] = df[col].astype(str).str.replace('%', '').astype(float)
            print(f"已将 {col} 转换为数值")
        except Exception as e:
            print(f"转换 {col} 失败: {e}")
            sys.exit(1)
    
    return df

def create_single_plot(df, output_prefix, width=14, height=8, dpi=300, color_palette='Set2'):
    """创建单个大图"""
    
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'Arial']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(width, height))
    
    # 创建分组条形图
    n_groups = len(df)
    bar_width = 0.2
    index = range(n_groups)
    
    # 设置颜色
    colors = plt.cm.get_cmap(color_palette, 4)
    
    # 绘制四个指标
    bar1 = ax.bar([i - 1.5*bar_width for i in index], 
                   df['long∩up'], bar_width, 
                   label='long∩up', color=colors(0), edgecolor='black')
    bar2 = ax.bar([i - 0.5*bar_width for i in index], 
                   df['long∩down'], bar_width, 
                   label='long∩down', color=colors(1), edgecolor='black')
    bar3 = ax.bar([i + 0.5*bar_width for i in index], 
                   df['short∩up'], bar_width, 
                   label='short∩up', color=colors(2), edgecolor='black')
    bar4 = ax.bar([i + 1.5*bar_width for i in index], 
                   df['short∩down'], bar_width, 
                   label='short∩down', color=colors(3), edgecolor='black')
    
    # 添加标签和标题
    ax.set_xlabel('compare', fontsize=12)
    ax.set_ylabel('ratio(%)', fontsize=12)
    ax.set_title('DUAPA vs DEAPA', fontsize=14, fontweight='bold')
    
    # 设置x轴标签
    if 'group1' in df.columns:
        # 使用group1和compare组合作为标签
        x_labels = [f"{row['group1']}\n{row['compare']}" for _, row in df.iterrows()]
    else:
        x_labels = df['compare']
    
    ax.set_xticks(index)
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=10)
    ax.legend(fontsize=10)
    
    # 添加网格
    ax.grid(True, axis='y', alpha=0.3, linestyle='--')
    
    # 调整y轴范围，留出空间
    y_max = max(df[['long∩up', 'long∩down', 'short∩up', 'short∩down']].max().max(), 10)
    ax.set_ylim(0, y_max * 1.1)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存为PDF
    pdf_file = f"{output_prefix}_single.pdf"
    plt.savefig(pdf_file, dpi=dpi, format='pdf', bbox_inches='tight')
    print(f"单个大图PDF已保存: {pdf_file}")
    
    # 保存为PNG
    png_file = f"{output_prefix}_single.png"
    plt.savefig(png_file, dpi=dpi, format='png', bbox_inches='tight')
    print(f"单个大图PNG已保存: {png_file}")
    
    plt.close(fig)

def create_grouped_plots(df, output_prefix, width=14, height=8, dpi=300, color_palette='Set2'):
    """按group1分组创建多个子图"""
    
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'Arial']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 获取所有group1
    groups = df['group1'].unique()
    n_groups = len(groups)
    
    # 计算子图布局
    cols = min(3, n_groups)
    rows = (n_groups + cols - 1) // cols
    
    # 创建图形
    fig, axes = plt.subplots(rows, cols, figsize=(width, height), squeeze=False)
    fig.suptitle('DUAPA vs DEAPA (groups)', fontsize=16, fontweight='bold')
    
    # 设置颜色
    colors = plt.cm.get_cmap(color_palette, 4)
    
    for idx, group in enumerate(groups):
        row = idx // cols
        col = idx % cols
        ax = axes[row, col]
        
        # 获取当前组的数据
        group_df = df[df['group1'] == group]
        
        # 创建当前组的图表
        n_bars = len(group_df)
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
        
        # 设置当前子图标题和标签
        ax.set_title(f'{group}', fontsize=12, fontweight='bold')
        ax.set_xlabel('compare', fontsize=10)
        ax.set_ylabel('ratio(%)', fontsize=10)
        
        # 设置x轴标签
        x_labels = group_df['compare']
        ax.set_xticks(index)
        ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
        
        # 添加网格
        ax.grid(True, axis='y', alpha=0.3, linestyle='--')
        
        # 调整y轴范围
        y_max = max(group_df[['long∩up', 'long∩down', 'short∩up', 'short∩down']].max().max(), 10)
        ax.set_ylim(0, y_max * 1.1)
        
        # 仅第一个子图显示图例
        if idx == 0:
            ax.legend(fontsize=8)
    
    # 隐藏多余的子图
    for idx in range(n_groups, rows * cols):
        row = idx // cols
        col = idx % cols
        axes[row, col].axis('off')
    
    # 调整布局
    plt.tight_layout()
    
    # 保存为PDF
    pdf_file = f"{output_prefix}_grouped.pdf"
    plt.savefig(pdf_file, dpi=dpi, format='pdf', bbox_inches='tight')
    print(f"分组PDF已保存: {pdf_file}")
    
    # 保存为PNG
    png_file = f"{output_prefix}_grouped.png"
    plt.savefig(png_file, dpi=dpi, format='png', bbox_inches='tight')
    print(f"分组PNG已保存: {png_file}")
    
    plt.close(fig)

def create_stacked_plot(df, output_prefix, width=14, height=8, dpi=300, color_palette='Set2'):
    """创建堆叠条形图"""
    
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'Arial']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 设置颜色
    colors = plt.cm.get_cmap(color_palette, 4)
    
    # 创建图形
    fig, axes = plt.subplots(1, 2, figsize=(width, height))
    fig.suptitle('DUAPA vs DEAPA - stacked', fontsize=16, fontweight='bold')
    
    # 准备数据
    categories = [f"{row['group1']}\n{row['compare']}" for _, row in df.iterrows()]
    
    # 左侧：long相关数据堆叠
    axes[0].bar(categories, df['long∩up'], label='long∩up', 
                color=colors(0), edgecolor='black')
    axes[0].bar(categories, df['long∩down'], bottom=df['long∩up'], 
                label='long∩down', color=colors(1), edgecolor='black')
    
    axes[0].set_title('Lengthen', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('ratio(%)', fontsize=10)
    axes[0].set_xticklabels(categories, rotation=45, ha='right', fontsize=8)
    axes[0].legend(fontsize=9)
    axes[0].grid(True, axis='y', alpha=0.3, linestyle='--')
    
    # 右侧：short相关数据堆叠
    axes[1].bar(categories, df['short∩up'], label='short∩up', 
                color=colors(2), edgecolor='black')
    axes[1].bar(categories, df['short∩down'], bottom=df['short∩up'], 
                label='short∩down', color=colors(3), edgecolor='black')
    
    axes[1].set_title('Shorten', fontsize=12, fontweight='bold')
    axes[1].set_ylabel('ratio(%)', fontsize=10)
    axes[1].set_xticklabels(categories, rotation=45, ha='right', fontsize=8)
    axes[1].legend(fontsize=9)
    axes[1].grid(True, axis='y', alpha=0.3, linestyle='--')
    
    # 调整y轴范围
    for ax in axes:
        y_max = ax.get_ylim()[1]
        ax.set_ylim(0, y_max * 1.1)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存
    pdf_file = f"{output_prefix}_stacked.pdf"
    plt.savefig(pdf_file, dpi=dpi, format='pdf', bbox_inches='tight')
    print(f"堆叠图PDF已保存: {pdf_file}")
    
    png_file = f"{output_prefix}_stacked.png"
    plt.savefig(png_file, dpi=dpi, format='png', bbox_inches='tight')
    print(f"堆叠图PNG已保存: {png_file}")
    
    plt.close(fig)

def create_heatmap(df, output_prefix, width=14, height=8, dpi=300):
    """创建热图"""
    
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'Arial']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 准备数据
    heatmap_data = df[['long∩up', 'long∩down', 'short∩up', 'short∩down']].values.T
    row_labels = ['long∩up', 'long∩down', 'short∩up', 'short∩down']
    col_labels = [f"{row['group1']} {row['compare']}" for _, row in df.iterrows()]
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(width, height))
    
    # 创建热图
    im = ax.imshow(heatmap_data, cmap='YlOrRd', aspect='auto')
    
    # 设置刻度
    ax.set_xticks(range(len(col_labels)))
    ax.set_yticks(range(len(row_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha='right', fontsize=9)
    ax.set_yticklabels(row_labels, fontsize=10)
    
    # 添加颜色条
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('ratio(%)', rotation=-90, va="bottom", fontsize=10)
    
    # 添加标题
    ax.set_title('DUAPA vs DEAPA Heatmap', fontsize=14, fontweight='bold')
    
    # 在单元格中添加数值
    for i in range(len(row_labels)):
        for j in range(len(col_labels)):
            text = ax.text(j, i, f'{heatmap_data[i, j]:.1f}',
                          ha="center", va="center", color="black", fontsize=8)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存
    pdf_file = f"{output_prefix}_heatmap.pdf"
    plt.savefig(pdf_file, dpi=dpi, format='pdf', bbox_inches='tight')
    print(f"热图PDF已保存: {pdf_file}")
    
    png_file = f"{output_prefix}_heatmap.png"
    plt.savefig(png_file, dpi=dpi, format='png', bbox_inches='tight')
    print(f"热图PNG已保存: {png_file}")
    
    plt.close(fig)

def main():
    """主函数"""
    args = parse_arguments()
    
    # 加载数据
    df = load_data(args.input_file)
    
    # 验证和转换数据
    df = validate_data(df)
    
    # 按group1排序
    if 'group1' in df.columns:
        df = df.sort_values(['group1', 'compare'])
    
    print(f"\n数据统计:")
    print(df[['long∩up', 'long∩down', 'short∩up', 'short∩down']].describe())
    
    # 创建输出目录（如果需要）
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 根据样式参数创建图表
    if args.style == 'single':
        print("\n正在创建单个图表...")
        create_single_plot(df, args.output_prefix, args.width, args.height, args.dpi, args.color_palette)
    elif args.style == 'grouped':
        print("\n正在创建分组图表...")
        create_grouped_plots(df, args.output_prefix, args.width, args.height, args.dpi, args.color_palette)
    elif args.style == 'stacked':
        print("\n正在创建堆叠图表...")
        create_stacked_plot(df, args.output_prefix, args.width, args.height, args.dpi, args.color_palette)
    elif args.style == 'heatmap':
        print("\n正在创建热图...")
        create_heatmap(df, args.output_prefix, args.width, args.height, args.dpi)
    elif args.style == 'all':
        print("\n正在创建所有图表...")
        create_single_plot(df, args.output_prefix, args.width, args.height, args.dpi, args.color_palette)
        create_grouped_plots(df, args.output_prefix, args.width, args.height, args.dpi, args.color_palette)
        create_stacked_plot(df, args.output_prefix, args.width, args.height, args.dpi, args.color_palette)
        create_heatmap(df, args.output_prefix, args.width, args.height, args.dpi)
    else:
        print(f"未知的样式: {args.style}")
        print("使用 --style single|grouped|stacked|heatmap|all")
        sys.exit(1)
    
    print(f"\n所有图表已保存到: {args.output_prefix}*")
    print("完成！")

if __name__ == "__main__":
    main()
