import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import argparse
import sys

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='生成基因表达热图')
    parser.add_argument('-e', '--expression', required=True, help='表达量矩阵文件(TSV)')
    parser.add_argument('-g', '--groups', required=True, help='样本分组文件(TSV)')
    parser.add_argument('-i', '--genes', help='基因ID列表文件(可选)')
    parser.add_argument('-o', '--output', default='gene_expression_heatmap.pdf', 
                       help='输出PDF文件名')
    parser.add_argument('--dpi', type=int, default=300, help='输出图像DPI')
    args = parser.parse_args()

    try:
        # 读取表达量矩阵
        expr_matrix = pd.read_csv(args.expression, sep='\t', index_col=0)
        
        # 检查表达矩阵是否为空
        if expr_matrix.empty:
            raise ValueError("表达矩阵为空")

        # 读取分组信息
        group_info = pd.read_csv(args.groups, sep='\t', header=None, 
                               names=['sample', 'group'])
        
        # 检查分组文件是否为空
        if group_info.empty:
            raise ValueError("分组文件为空")

        # 获取两组共有的样本
        common_samples = list(set(group_info['sample']) & set(expr_matrix.columns))
        if not common_samples:
            raise ValueError("表达矩阵和分组文件没有共同的样本ID")
        
        # 筛选样本和排序
        group_info = group_info[group_info['sample'].isin(common_samples)]
        
        # 按分组排序样本
        group_info = group_info.sort_values('group')
        common_samples = group_info['sample'].tolist()
        expr_matrix = expr_matrix[common_samples]

        # 如果有提供基因列表文件，筛选基因
        if args.genes:
            with open(args.genes) as f:
                target_genes = [line.strip() for line in f if line.strip()]
            expr_matrix = expr_matrix.loc[expr_matrix.index.intersection(target_genes)]
            
            if expr_matrix.empty:
                raise ValueError("没有找到基因列表中的基因ID")

        # 创建颜色映射
        unique_groups = group_info['group'].unique()
        palette = sns.color_palette("husl", len(unique_groups))
        group_colors = dict(zip(unique_groups, palette))

        # 创建列颜色条
        col_colors = group_info.set_index('sample')['group'].map(group_colors).loc[common_samples]

        # 计算图形高度(基于基因数量)
        fig_height = max(6, len(expr_matrix) * 0.2)
        fig_width = max(8, len(common_samples) * 0.1)

        # 创建图形
        fig, (cbar_ax, heatmap_ax) = plt.subplots(
            2, 1, gridspec_kw={"height_ratios": [0.02, 1]}, 
            figsize=(fig_width, fig_height))
        
        # 使用对数归一化
        norm = Normalize(vmin=np.percentile(expr_matrix.values, 5),
                        vmax=np.percentile(expr_matrix.values, 95))

        # 绘制热图
        heatmap = sns.heatmap(
            expr_matrix,
            cmap='viridis',
            norm=norm,
            yticklabels=True,
            xticklabels=False,
            ax=heatmap_ax,
            cbar_ax=cbar_ax,
            cbar_kws={"orientation": "horizontal"})

        # 添加分组颜色条
        for i, (sample, group) in enumerate(zip(common_samples, group_info['group'])):
            heatmap_ax.add_patch(plt.Rectangle(
                (i, -0.05), 1, 0.05, 
                color=group_colors[group],
                transform=heatmap_ax.get_xaxis_transform(),
                clip_on=False))

        # 添加分组图例
        handles = [plt.Rectangle((0,0),1,1, color=group_colors[g]) for g in unique_groups]
        heatmap_ax.legend(
            handles, unique_groups,
            title='Group',
            bbox_to_anchor=(1.02, 1),
            loc='upper left',
            borderaxespad=0.)

        # 调整基因名称字体大小
        heatmap_ax.set_yticklabels(
            heatmap_ax.get_yticklabels(),
            fontsize=6)

        # 调整布局
        plt.tight_layout()

        # 保存为PDF
        plt.savefig(args.output, bbox_inches='tight', dpi=args.dpi)
        print(f"热图已保存为 {args.output}")
        
    except Exception as e:
        print(f"错误发生: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
