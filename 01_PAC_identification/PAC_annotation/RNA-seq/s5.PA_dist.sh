pro=$1
bin=//MZAPA/ZY/post_analysis/PAC_CS_RNAseq_feature/bin
dir=//MZAPA/ZY/post_analysis/PAC_CS_RNAseq_feature
cd $dir/$pro/filter01

#计算pA的位置分布
grep three_prime_UTR $pro\.intersect.allType|cut -f4,7-9,12|sort -u >tmp
perl $bin/distance_ratio.pl tmp $pro\.PA.3UTR_dis.txt
rm tmp
/share/apps/R-4.2.1/bin/Rscript $bin/density_plot_forRatio.r $pro\.PA.3UTR_dis.txt $pro\.PA.3UTR_dis

#计算PAC长度分布
grep three_prime_UTR $pro\.intersect.allType|awk '{print $5"\t"$3-$2+1}' |sort -u >$pro\.PAC.3UTR.len
sed -i '1i id\tlength' $pro\.PAC.3UTR.len
/share/apps/R-4.2.1/bin/Rscript $bin/length_bin_single.r $pro\.PAC.3UTR.len $pro\.PAC.3UTR.len

