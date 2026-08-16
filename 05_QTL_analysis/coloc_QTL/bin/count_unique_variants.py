#!/usr/bin/env python3
"""
统计输入文件指定列中逗号分隔的位点唯一总数（去重后）。
默认列索引为5（从1开始计数），默认跳过第一行表头。
"""

import argparse
import sys
from collections import Counter

def main():
    parser = argparse.ArgumentParser(description='统计指定列中逗号分隔位点的唯一总数')
    parser.add_argument('-i', '--input', required=True, help='输入文件 (制表符分隔)')
    parser.add_argument('-o', '--output', default=None,
                        help='可选：输出每个位点及其出现次数的文件 (制表符)')
    parser.add_argument('--col', type=int, default=5,
                        help='目标列索引 (从1开始计数, 默认5)')
    parser.add_argument('--skip_header', action='store_true', default=True,
                        help='跳过第一行 (默认开启)')
    parser.add_argument('--no_skip_header', dest='skip_header', action='store_false',
                        help='不跳过第一行')
    parser.add_argument('--delimiter', default='\t',
                        help='列分隔符 (默认制表符)')
    args = parser.parse_args()

    unique_variants = set()

    with open(args.input, 'r') as f:
        lines = f.readlines()
        if args.skip_header and lines:
            lines = lines[1:]  # 跳过表头

        for line in lines:
            line = line.strip()
            if not line:
                continue
            cols = line.split(args.delimiter)
            if len(cols) < args.col:
                # 该行没有足够的列，跳过
                continue
            # 提取目标列 (索引从0开始)
            cell = cols[args.col - 1].strip()
            if not cell:
                continue
            # 按逗号分割
            variants = [v.strip() for v in cell.split(',') if v.strip()]
            unique_variants.update(variants)

    total = len(unique_variants)
    print(f"总唯一位点数: {total}")

    if args.output:
        # 重新读取文件统计每个位点的出现次数
        counter = Counter()
        with open(args.input, 'r') as f:
            lines = f.readlines()
            if args.skip_header and lines:
                lines = lines[1:]
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                cols = line.split(args.delimiter)
                if len(cols) < args.col:
                    continue
                cell = cols[args.col - 1].strip()
                if not cell:
                    continue
                variants = [v.strip() for v in cell.split(',') if v.strip()]
                counter.update(variants)
        with open(args.output, 'w') as out:
            for variant, count in sorted(counter.items()):
                out.write(f"{variant}\t{count}\n")
        print(f"每个位点的计数已保存到: {args.output}")

if __name__ == "__main__":
    main()
