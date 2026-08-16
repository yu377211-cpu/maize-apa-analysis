for i in ln mixed nss popcorn ss sweet ts
do
	bash ../bin/updown_elongeshort_list.sh ${i} >${i}.updown_longshort.list
done

for i in ln mixed nss popcorn ss sweet ts
do
	bash ../bin/updown_elongeshort.sh ${i} >${i}.updown_longshort.num.txt
done

Rscript ../bin/combine_matrix.r
统计信息:
总样本数: 126
总特征数: 1161
标记为1的数量: 2599
标记为-1的数量: 3176
标记为0的数量: 140511

文件已保存:
  feature_matrix.csv - 原始矩阵文件
  feature_matrix_with_groups.csv - 包含分组信息
  feature_matrix_transposed.csv - 转置矩阵
✔ 成功生成 main_matrix_by_subpop.csv - 聚合矩阵文件
✔ 成功生成 supp_matrix_by_subpop_shared.csv - 比较组的矩阵文件

Rscript ../bin/heatmap_matrix.r

cat *.list >Total.updown_longshort_list  #并删除多余的行，添加表头
python3 ../bin/calculate.final_score.py -i Total.updown_longshort_list  -o Total.updown_elongeshort.score
Rscript ../bin/plot_final_score_heatmap.R -i Total.updown_elongeshort.score -o APA_coupling_subpop -t 30 -m 2
Rscript ../bin/coupling_consistency_dotplot.R -i Total.updown_longshort_list -o same_subpop.updown_longshort
python3 ../bin/calc_asymmetry_DCBindex.py -i "*.updown_longshort.list" -o subpop_asymmetry

python3 ../bin/merge_updown.py >sameSubpop.NUMsta.txt
wc -l /APA/AFreview/PAC_3seq_diff_usage/same_subpop/*/*.sig.lengthening >lengthening
wc -l /APA/AFreview/PAC_3seq_diff_usage/same_subpop/*/*.sig.shortening >shortening
sed -i 's/^[[:space:]]*//' lengthening
sed -i 's/^[[:space:]]*//' shortening
==> sameSubpop.lengthening_shortening.sta.xls
==> awk -F'\t' 'NR==1 || $3!="0.00%" || $4!="0.00%" || $5!="0.00%" || $6!="0.00%"' sameSubpop.ratio.txt
