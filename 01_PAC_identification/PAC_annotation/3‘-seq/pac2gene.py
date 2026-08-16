#!/usr/bin/env python3
"""
将 PAC 水平的 read count 矩阵合并为基因水平的 read count 矩阵（求和）。
用法：
    python pac_to_gene.py -p pac_counts.tsv -g gene_pac_map.tsv -o gene_counts.tsv
"""

import argparse
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(description="Merge PAC-level read counts to gene-level by summing counts per gene.")
    parser.add_argument("-p", "--pac_matrix", required=True, help="PAC count matrix (TSV), first column = PAC ID, columns = samples")
    parser.add_argument("-g", "--gene_map", required=True, help="Gene-PAC mapping file (TSV), first column = gene ID, other columns = PAC IDs")
    parser.add_argument("-o", "--output", required=True, help="Output gene count matrix (TSV)")
    args = parser.parse_args()

    # 1. 读取基因-PAC 映射
    print("Reading gene-PAC mapping...", file=sys.stderr)
    gene_map = {}
    with open(args.gene_map, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            gene = parts[0]
            pacs = [p for p in parts[1:] if p]  # 去除可能的空值
            gene_map[gene] = pacs
    print(f"Loaded {len(gene_map)} genes.", file=sys.stderr)

    # 2. 读取 PAC 计数矩阵（第一列为 PAC ID，作为行索引）
    print("Reading PAC count matrix...", file=sys.stderr)
    pac_df = pd.read_csv(args.pac_matrix, sep='\t', index_col=0)
    print(f"PAC matrix shape: {pac_df.shape}", file=sys.stderr)

    # 检查 PAC ID 是否在矩阵中，过滤掉缺失的
    all_pacs = set(pac_df.index)
    gene_results = {}
    for gene, pac_list in gene_map.items():
        # 取在矩阵中存在的 PAC
        valid_pacs = [p for p in pac_list if p in all_pacs]
        if not valid_pacs:
            # 若无任何 PAC 在矩阵中，可输出全零行（或跳过）
            # 此处跳过（不输出该基因）
            continue
        # 提取这些 PAC 的行，按列求和
        gene_counts = pac_df.loc[valid_pacs].sum(axis=0)
        gene_results[gene] = gene_counts

    # 3. 构建结果 DataFrame
    if not gene_results:
        print("No genes with valid PACs found. Output empty.", file=sys.stderr)
        # 创建一个空 DataFrame 以保证输出格式
        result_df = pd.DataFrame(columns=pac_df.columns)
    else:
        result_df = pd.DataFrame.from_dict(gene_results, orient='index')
        # 确保列顺序与 PAC 矩阵一致
        result_df = result_df[pac_df.columns]
        # 基因名为行名，第一列需要输出为基因 ID
        result_df.index.name = 'Gene'

    # 4. 写入输出
    print(f"Writing gene count matrix ({result_df.shape[0]} genes, {result_df.shape[1]} samples)...", file=sys.stderr)
    result_df.to_csv(args.output, sep='\t', header=True, index=True)
    print("Done.", file=sys.stderr)

if __name__ == "__main__":
    main()
