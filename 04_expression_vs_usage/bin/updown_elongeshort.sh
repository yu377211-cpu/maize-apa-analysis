#!/bin/bash
group="$1"
#DIFF_DIR="/share/pub/xingsl/shilai/project/APA/AFreview/PAC_3seq_diffana/PAC_readcount/same_subpop/diff/$group"
#APA_DIR="/share/pub/xingsl/shilai/project/APA/AFreview/PAC_3seq_diff_usage/same_subpop/$group"
DIFF_DIR="/share/pub/xingsl/shilai/project/APA/AFreview/PAC_3seq_diffana/PAC_readcount/same_tissue/diff/$group"
APA_DIR="/share/pub/xingsl/shilai/project/APA/AFreview/PAC_3seq_diff_usage/same_tissue/$group"

echo "比较组              上调∩long  下调∩long   上调∩short   下调∩short"
echo "------------------------------------------------------------------"

cd "$DIFF_DIR" 2>/dev/null || { echo "无法进入目录: $DIFF_DIR"; exit 1; }

# 安全的 grep -c：无论什么情况都只输出一个干净的 0
safe_count() {
    local n
    n=$(grep -c -f "$1" "$2" 2>/dev/null)
    printf '%d' "${n:-0}"
}

for up_file in *.condition.diff.addgene.up; do
    comp="${up_file%.condition.diff.addgene.up}"
    up_path="$DIFF_DIR/${comp}.condition.diff.addgene.up"
    down_path="$DIFF_DIR/${comp}.condition.diff.addgene.down"
    short_path="$APA_DIR/${group}.${comp}.sig.shortening"
    long_path="$APA_DIR/${group}.${comp}.sig.lengthening"

    # 关键的两个文件一定要存在；短/长文件可以缺/空
    if [[ -f "$short_path" && -f "$long_path" ]]; then
        c1=$(safe_count "$up_path"   "$long_path")
        c2=$(safe_count "$down_path" "$long_path")
        c3=$(safe_count "$up_path"   "$short_path")
        c4=$(safe_count "$down_path" "$short_path")
        printf "%-20s %10d %10d %10d %10d\n" "$comp" "$c1" "$c2" "$c3" "$c4"
    fi
done

echo "------------------------------------------------------------------"

