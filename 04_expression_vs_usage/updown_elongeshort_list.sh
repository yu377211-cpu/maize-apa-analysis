#!/bin/bash

group="LMAN"
DIFF_DIR="/APA/PAC_3seq_diffana/PAC_tpm_ab22/same_tissue/diff/$group"
APA_DIR="/APA/PAC_3seq_diff_usage/same_tissue/$group"

echo "组   比较组              特征ID                类型   标记"
echo "--------------------------------------------------------------"

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
    
    # 检查APA差异文件是否存在
    if [[ ! -f "$short_path" || ! -f "$long_path" ]]; then
        printf "%-5s %-20s %-20s %-10s %5s\n" "$group" "$comp" "" "" "NA"
        continue
    fi
    
    # 检查差异表达文件是否存在
    if [[ ! -f "$up_path" || ! -f "$down_path" ]]; then
        printf "%-5s %-20s %-20s %-10s %5s\n" "$group" "$comp" "" "" "NA"
        continue
    fi
    
    # 记录是否有重叠特征
    has_overlap=false
    
    # 处理long-up类型 (标记1)
    if [[ -f "$up_path" && -f "$long_path" ]]; then
        grep -f "$up_path" "$long_path" 2>/dev/null | while read -r feature; do
            if [[ -n "$feature" ]]; then
                printf "%-5s %-20s %-20s %-10s %5d\n" "$group" "$comp" "$feature" "long-up" 1
                has_overlap=true
            fi
        done
    fi
    
    # 处理short-down类型 (标记1)
    if [[ -f "$down_path" && -f "$short_path" ]]; then
        grep -f "$down_path" "$short_path" 2>/dev/null | while read -r feature; do
            if [[ -n "$feature" ]]; then
                printf "%-5s %-20s %-20s %-10s %5d\n" "$group" "$comp" "$feature" "short-down" 1
                has_overlap=true
            fi
        done
    fi
    
    # 处理long-down类型 (标记-1)
    if [[ -f "$down_path" && -f "$long_path" ]]; then
        grep -f "$down_path" "$long_path" 2>/dev/null | while read -r feature; do
            if [[ -n "$feature" ]]; then
                printf "%-5s %-20s %-20s %-10s %5d\n" "$group" "$comp" "$feature" "long-down" -1
                has_overlap=true
            fi
        done
    fi
    
    # 处理short-up类型 (标记-1)
    if [[ -f "$up_path" && -f "$short_path" ]]; then
        grep -f "$up_path" "$short_path" 2>/dev/null | while read -r feature; do
            if [[ -n "$feature" ]]; then
                printf "%-5s %-20s %-20s %-10s %5d\n" "$group" "$comp" "$feature" "short-up" -1
                has_overlap=true
            fi
        done
    fi
    
    # 如果没有找到任何重叠的特征，但文件存在
    # 注意：这里需要等while循环结束后再判断
    # 由于while循环在子进程中运行，我们需要用其他方法跟踪是否找到重叠
    
    # 创建一个临时方法来检查是否有重叠
    check_overlap() {
        local has_overlap=false
        
        # 检查long-up
        if [[ -f "$up_path" && -f "$long_path" ]]; then
            if grep -q -f "$up_path" "$long_path" 2>/dev/null; then
                has_overlap=true
            fi
        fi
        
        # 检查short-down
        if [[ -f "$down_path" && -f "$short_path" ]]; then
            if grep -q -f "$down_path" "$short_path" 2>/dev/null; then
                has_overlap=true
            fi
        fi
        
        # 检查long-down
        if [[ -f "$down_path" && -f "$long_path" ]]; then
            if grep -q -f "$down_path" "$long_path" 2>/dev/null; then
                has_overlap=true
            fi
        fi
        
        # 检查short-up
        if [[ -f "$up_path" && -f "$short_path" ]]; then
            if grep -q -f "$up_path" "$short_path" 2>/dev/null; then
                has_overlap=true
            fi
        fi
        
        echo "$has_overlap"
    }
    
    # 检查是否有重叠
    overlap_result=$(check_overlap)
    
    if [[ "$overlap_result" == "false" ]]; then
        printf "%-5s %-20s %-20s %-10s %5d\n" "$group" "$comp" "" "" 0
    fi
done

echo "--------------------------------------------------------------"
