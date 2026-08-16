#!/bin/bash

group="LMAN"
DIFF_DIR="/APA/PAC_3seq_diffana/PAC_tpm_ab22/same_tissue/diff/$group"
APA_DIR="/APA/PAC_3seq_diff_usage/same_tissue/$group"

echo "比较组              上调∩long  下调∩long   上调∩short   下调∩short"
echo "------------------------------------------------------------------"

# 直接进入目录处理
cd "$DIFF_DIR" 2>/dev/null || { echo "无法进入目录: $DIFF_DIR"; exit 1; }

for up_file in *.condition.TPM.diff_sig.up; do
    # 提取比较名称
    comp="${up_file%.condition.TPM.diff_sig.up}"
    
    # 构建文件路径
    up_path="$DIFF_DIR/${comp}.condition.TPM.diff_sig.up"
    down_path="$DIFF_DIR/${comp}.condition.TPM.diff_sig.down"
    short_path="$APA_DIR/${group}.${comp}.sig.shortening"
    long_path="$APA_DIR/${group}.${comp}.sig.lengthening"
    
    # 检查文件是否存在
    if [[ -f "$short_path" && -f "$long_path" ]]; then
        # 执行统计
        c3=$(grep -c -f "$up_path" "$short_path" 2>/dev/null || echo 0)
        c1=$(grep -c -f "$up_path" "$long_path" 2>/dev/null || echo 0)
        c2=$(grep -c -f "$down_path" "$long_path" 2>/dev/null || echo 0)
        c4=$(grep -c -f "$down_path" "$short_path" 2>/dev/null || echo 0)
        
        printf "%-20s %10d %10d %10d %10d\n" "$comp" "$c1" "$c2" "$c3" "$c4"
    fi
done

echo "------------------------------------------------------------------"
