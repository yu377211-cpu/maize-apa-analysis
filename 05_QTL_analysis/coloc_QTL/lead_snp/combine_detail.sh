# 获取所有文件
files=(*_coloc.detail.txt)
# 先输出表头
echo -e "Tissue\tPAC\tae_cs\ttranscript\ta_cs\tshared_variants" > combine.coloc.detail.txt
# 对每个文件，跳过第一行，打印组织名和内容
for f in "${files[@]}"; do
    tissue=$(basename "$f" _coloc.detail.txt)
    awk -v t="$tissue" 'NR>1 {print t "\t" $0}' "$f"
done >> combine.coloc.detail.txt
