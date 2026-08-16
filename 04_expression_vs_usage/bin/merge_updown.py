#!/usr/bin/env python3
import glob
import re
import sys

def process_file(filepath):
    """解析单个文件，返回 (group, rows)"""
    group = filepath.split('.')[0]          # 取第一个点之前作为group1
    rows = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # 跳过表头（包含“比较组”或“上调”等）
            if '比较组' in line or '上调∩long' in line:
                continue
            # 跳过分隔线：全是 '-' 或 '--' 开头（去除空格后）
            stripped = re.sub(r'\s+', '', line)   # 去掉所有空白
            if stripped.startswith('--') or stripped.startswith('---'):
                continue
            # 按空白符（空格/制表符）切分
            parts = re.split(r'\s+', line)
            # 过滤空字符串（防止分割出空项）
            parts = [p for p in parts if p != '']
            if len(parts) < 5:
                continue   # 列数不足，跳过
            compare_raw = parts[0]
            # 将 'vs' 替换为 ' vs ' （前后空格）
            compare_new = compare_raw.replace('vs', ' vs ')
            # 后四列
            values = parts[1:5]
            rows.append((compare_new, values))
    return group, rows

def main():
    files = glob.glob('*.updown_longshort.num.txt')
    if not files:
        print("未找到任何 .updown_longshort.num.txt 文件", file=sys.stderr)
        sys.exit(1)

    # 输出表头（制表符分隔）
    header = ['group1', 'compare', 'up∩long', 'down∩long', 'up∩short', 'down∩short']
    print('\t'.join(header))

    for f in sorted(files):
        group, rows = process_file(f)
        for compare_str, vals in rows:
            out_line = [group, compare_str] + vals
            print('\t'.join(out_line))

if __name__ == '__main__':
    main()
