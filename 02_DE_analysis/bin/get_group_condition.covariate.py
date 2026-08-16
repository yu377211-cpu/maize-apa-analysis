#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import re
import pandas as pd
import argparse

__author__ = 'zhangyu'
__date__ = '2021-08-05'

def get_group(compare_name):
    """从比较名称拆分出两个组名"""
    group_tre, group_ck = re.split('vs', compare_name)
    return group_tre, group_ck

def set_list(list_a, list_b):
    """返回交集列表"""
    return list(set(list_a).intersection(set(list_b)))

def judge_compare(pd_condition, group_a, group_b):
    """
    在数据框中找到包含 group_a 和 group_b 的列（通常为同一列），
    筛选该列值等于 group_a 或 group_b 的所有行，返回完整 DataFrame 的子集。
    """
    # 查找包含 group_a 的列
    loc_a = (pd_condition == group_a).any()
    if not loc_a.any():
        print(f"ERROR: group '{group_a}' not found in any column.")
        sys.exit(1)
    col_a = loc_a.index[loc_a][0]

    # 查找包含 group_b 的列
    loc_b = (pd_condition == group_b).any()
    if not loc_b.any():
        print(f"ERROR: group '{group_b}' not found in any column.")
        sys.exit(1)
    col_b = loc_b.index[loc_b][0]

    if col_a == col_b:
        # 同一列，直接筛选
        mask = (pd_condition[col_a] == group_a) | (pd_condition[col_b] == group_b)
        filtered_df = pd_condition[mask].copy()
    else:
        # 不同列，分别筛选后合并，并检查样本重复
        df_a = pd_condition[pd_condition[col_a] == group_a]
        df_b = pd_condition[pd_condition[col_b] == group_b]
        duplicate_idx = df_a.index.intersection(df_b.index)
        if not duplicate_idx.empty:
            print(f"ERROR: group_a '{group_a}' and group_b '{group_b}' have duplicate samples: {duplicate_idx.tolist()}")
            sys.exit(1)
        filtered_df = pd.concat([df_a, df_b])

    return filtered_df

def write_to_file(df, outfile):
    """将 DataFrame 写入 TSV 文件，自动创建目录，保留索引（sample）和所有列"""
    outdir = os.path.dirname(os.path.realpath(outfile))
    if outdir and not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except FileExistsError:
            print(f"dir {outdir} already exists, not duplicated creating")
        except Exception:
            os.system(f"mkdir -p {outdir}")
    # 写入时包含索引，并将索引列命名为 'sample'
    df.to_csv(outfile, sep='\t', mode='w', index=True, index_label='sample')

def main():
    parser = argparse.ArgumentParser(
        description='Extract rows for two groups from a group file based on a comparison name.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f'author:\t{__author__}\ndate:\t{__date__}\n'
    )

    parser.add_argument('-c', '--compare', help='compare name, e.g. GRootvsGShoot', dest='compare', required=True)
    parser.add_argument('-o', '--outfile', help='output file path', metavar='', dest='outfile', required=True)
    parser.add_argument('-i', '--infile', help='input group file (TSV, first column as sample index)', metavar='',
                        dest='infile', required=True)

    args = parser.parse_args()

    # 读取文件，第一列为索引（sample）
    complex_condition = pd.read_csv(args.infile, header=0, sep='\t', index_col=0)
    group_tre, group_ck = get_group(args.compare)
    filtered_df = judge_compare(complex_condition, group_tre, group_ck)
    write_to_file(filtered_df, args.outfile)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    main()
