#!/bin/bash

# 检查是否提供了输入文件作为参数
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

# 读取文件名
filename=$1

# 检查文件是否存在
if [ ! -f "$filename" ]; then
    echo "File not found: $filename"
    exit 1
fi

# 使用awk处理文件
awk 'BEGIN {
    FS = OFS = "\t"
}

NR == 1 {
    # 读取标题
    for (i = 9; i <= NF; i++) {
        titles[i - 8] = $i
    }
    next
}

NR > 1 {
    # 统计非零元素
    for (i = 9; i <= NF; i++) {
        if ($i != 0) {
            non_zero[titles[i - 8]]++
        }
    }
}

END {
    # 输出结果
    for (title in titles) {
        print titles[title], non_zero[titles[title]]
    }
}' "$filename" > $filename\.sampleSta
