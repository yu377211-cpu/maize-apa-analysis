import pandas as pd
import argparse

# 设置命令行参数解析
parser = argparse.ArgumentParser(description='Annotate positions with nearest gene.')
parser.add_argument('positions_file', type=str, help='Input positions file in CSV format')
parser.add_argument('bed_file', type=str, help='Input BED file')
parser.add_argument('output_file', type=str, help='Output file for annotated positions')
args = parser.parse_args()
## python3 site_annotation.py all_tissue_trans.snp.bed /share2/pub/xingsl/xingsl/data_download/MZAPA/Ref/B73_V4/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gene.bed all_tissue_trans.snp.annot_gene.txt

# 读取位点文件和BED文件
#positions_df = pd.read_csv(args.positions_file, sep='\t', dtype={'Chr': str, 'pos': int}) #all_tissue_cis.snp.pos
positions_df = pd.read_csv(args.positions_file, sep='\t', header=None, names=['Chr', 'pos'], dtype={'Chr': str, 'pos': int}) #给没有tittle的列赋名，并读取为字符串
genes_bed_df = pd.read_csv(args.bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'gene_name', 'score', 'strand']) #Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gene.bed

# 添加一列用于存储找到的基因和距离
positions_df['nearest_gene'] = ''
positions_df['distance'] = -1
positions_df['direction'] = ''

for index, position in positions_df.iterrows():
    # 筛选相同染色体上的基因
    same_chr_genes = genes_bed_df[genes_bed_df['chr'] == position['Chr']].copy()
 
    # 找到位点左边最近的基因
    same_chr_genes.loc[:, 'end'] = pd.to_numeric(same_chr_genes['end'], errors='coerce')
    left_genes = same_chr_genes[same_chr_genes['end'] < position['pos']]
    if not left_genes.empty:
      min_idx = left_genes['end'].sub(position['pos']).abs().idxmin()
      if min_idx in left_genes.index:
        closest_left_gene = left_genes.loc[min_idx]
        left_gene_name = closest_left_gene['gene_name']
        left_distance = position['pos'] - closest_left_gene['end']
      else:
        left_gene_name = 'None'
        left_distance = float('inf')
    else:
      left_gene_name = 'None'
      left_distance = float('inf')
    
    # 找到位点右边最近的基因
    same_chr_genes.loc[:, 'start'] = pd.to_numeric(same_chr_genes['start'], errors='coerce')
    right_genes = same_chr_genes[same_chr_genes['start'] > position['pos']]
    if not right_genes.empty:
      min_idx = right_genes['start'].sub(position['pos']).abs().idxmin()
      if min_idx in right_genes.index:
        closest_right_gene = right_genes.loc[min_idx]
        right_gene_name = closest_right_gene['gene_name']
        right_distance = closest_right_gene['start'] - position['pos']
      else:
        right_gene_name = 'None'
        right_distance = float('inf')
    else:
      right_gene_name = 'None'
      right_distance = float('inf')
    
    # 确定最近的基因以及距离
    if left_distance < right_distance:
        positions_df.at[index, 'nearest_gene'] = left_gene_name
        positions_df.at[index, 'distance'] = left_distance
        positions_df.at[index, 'direction'] = 'upstream'
    else:
        positions_df.at[index, 'nearest_gene'] = right_gene_name
        positions_df.at[index, 'distance'] = right_distance
        positions_df.at[index, 'direction'] = 'downstream'
    # 如果位点位于基因区域中间
    between_genes = same_chr_genes[(same_chr_genes['start'] <= position['pos']) & (same_chr_genes['end'] >= position['pos'])]
    if not between_genes.empty:
        between_gene = between_genes.iloc[0]
        positions_df.at[index, 'nearest_gene'] = between_gene['gene_name']
        positions_df.at[index, 'distance'] = 0
        positions_df.at[index, 'direction'] = 'within'

# 输出结果至新的CSV文件
positions_df.to_csv(args.output_file, sep='\t', index=False)
