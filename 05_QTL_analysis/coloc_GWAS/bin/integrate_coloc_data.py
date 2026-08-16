#!/usr/bin/env python3
"""
合并 LD 矩阵、cis 关联文件和基因注释文件，生成包含 GWAS 和 QTL 信息的表格。
用法：
  python integrate_coloc_data.py --tissue GRoot --qtl_type 3aQTL \
      --ld output_500kb_ld_cis_significant.txt \
      --cis output_500kb.cis.tsv \
      --anno SRP115041_ab22.intersect.3UTR.annoGene \
      --output integrated_result.txt
"""

import pandas as pd
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='整合 LD、cis 和注释数据')
    parser.add_argument('--tissue', required=True, help='组织名称（如 GRoot）')
    parser.add_argument('--qtl_type', required=True, help='QTL 类型（如 3aQTL）')
    parser.add_argument('--ld', required=True, help='LD 矩阵文件（SNP_A, SNP_B, R2）')
    parser.add_argument('--cis', required=True, help='cis 关联文件（含 chr.x, pos.x, locus_id.x, variant_id.x, pip.x, cs_purity.x, pos.y, trait.y, distance, cs.x, cs_size.x）')
    parser.add_argument('--anno', required=True, help='基因注释文件（PAC, Gene）')
    parser.add_argument('--output', required=True, help='输出文件路径（制表符分隔）')
    args = parser.parse_args()

    # 读取文件
    try:
        ld_df = pd.read_csv(args.ld, sep='\t')
        cis_df = pd.read_csv(args.cis, sep='\t')
        anno_df = pd.read_csv(args.anno, sep='\t', header=None, names=['PAC', 'Gene'])
    except Exception as e:
        print(f"读取输入文件失败: {e}", file=sys.stderr)
        sys.exit(1)

    # 检查必需列
    required_ld = ['SNP_A', 'SNP_B', 'R2']
    required_cis = ['chr.x', 'pos.x', 'locus_id.x', 'variant_id.x', 'pip.x', 'cs_purity.x',
                    'pos.y', 'trait.y', 'distance', 'cs.x', 'cs_size.x']
    if not all(col in ld_df.columns for col in required_ld):
        print("LD 文件缺少必需列", file=sys.stderr)
        sys.exit(1)
    if not all(col in cis_df.columns for col in required_cis):
        print("cis 文件缺少必需列", file=sys.stderr)
        sys.exit(1)

    # 构建 SNP 标识（chr_pos）
    cis_df['SNP_A'] = cis_df['chr.x'] + '_' + cis_df['pos.x'].astype(str)
    cis_df['SNP_B'] = cis_df['chr.x'] + '_' + cis_df['pos.y'].astype(str)

    # 合并 LD 与 cis：根据 SNP_A 和 SNP_B 匹配
    merged = pd.merge(ld_df, cis_df, on=['SNP_A', 'SNP_B'], how='inner')

    if merged.empty:
        print("警告: 没有匹配到任何行", file=sys.stderr)
        # 创建空输出（含所有列）
        out_cols = ['Tissue', 'QTL_type', 'QTL_phenotype', 'QTL_Gene', 'QTL_variant',
                    'QTL_PIP', 'QTL_CS', 'CS_size', 'QTL_purity',
                    'GWAS_trait', 'GWAS_variant_B73V4', 'GWAS_variant_B73V3',
                    'GWAS_P', 'GWAS_Gene', 'LD_R2', 'Distance (bp)']
        pd.DataFrame(columns=out_cols).to_csv(args.output, sep='\t', index=False)
        return

    # 合并基因注释：通过 locus_id.x 匹配 anno 第一列 (PAC)
    merged = pd.merge(merged, anno_df, left_on='locus_id.x', right_on='PAC', how='left')
    merged['Gene'] = merged['Gene'].fillna('NA')

    # 构造输出列（按指定顺序）
    result = pd.DataFrame({
        'Tissue': args.tissue,
        'QTL_type': args.qtl_type,
        'QTL_phenotype': merged['locus_id.x'],
        'QTL_Gene': merged['Gene'],
        'QTL_variant': merged['variant_id.x'],
        'QTL_PIP': merged['pip.x'],
        'QTL_CS': merged['cs.x'],
        'CS_size': merged['cs_size.x'],
        'QTL_purity': merged['cs_purity.x'],
        'GWAS_trait': merged['trait.y'],
        'GWAS_variant_B73V4': merged['SNP_B'],
        'GWAS_variant_B73V3': 'NA',
        'GWAS_P': 'NA',
        'GWAS_Gene': 'NA',
        'LD_R2': merged['R2'],
        'Distance (bp)': merged['distance']
    })

    # 写入输出
    result.to_csv(args.output, sep='\t', index=False)
    print(f"输出成功: {args.output}, 共 {len(result)} 行")

if __name__ == '__main__':
    main()
