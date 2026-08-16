pro=${pro}
bin=/share/pub/xingsl/shilai/project/APA/PAC_3seq_hq/motif/bin
dir=/MZAPA/ZY/post_analysis/PAC_CS_RNAseq_feature/$pro/filter01
bedtools=/share/apps/bedtools2/bin/bedtools
mkdir $dir/DNA_motif
cd $dir/DNA_motif
#碱基分布
awk -F"\t" '{if ($7 != "."){OFS="\t";print}}' $dir/$pro\.intersect.3UTR|cut -f1-5,12|sort -u > $pro\.intersect.3UTR.bed
perl $bin/PA.up_down_strand.pl $pro\.intersect.3UTR.bed $pro\.3UTR.ud.bed
grep "+" $pro\.3UTR.ud.bed|grep upstream |cut -f1-3|sort -u >$pro\.3UTR.up_plus.bed
grep "+" $pro\.3UTR.ud.bed|grep downstream|cut -f1-3|sort -u >$pro\.3UTR.down_plus.bed
grep "-" $pro\.3UTR.ud.bed|grep upstream |cut -f1-3|sort -u >$pro\.3UTR.up_minus.bed
grep "-" $pro\.3UTR.ud.bed|grep downstream|cut -f1-3|sort -u >$pro\.3UTR.down_minus.bed

$bedtools getfasta -fi /MZAPA/Ref/B73_V4/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.3UTR.up_plus.bed >$pro\.3UTR.up_plus.fa
$bedtools getfasta -fi /MZAPA/Ref/B73_V4/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.3UTR.down_plus.bed >$pro\.3UTR.down_plus.fa
$bedtools getfasta -fi /MZAPA/Ref/B73_V4/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.3UTR.up_minus.bed >$pro\.3UTR.up_minus.fa
$bedtools getfasta -fi /MZAPA/Ref/B73_V4/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.3UTR.down_minus.bed >$pro\.3UTR.down_minus.fa
perl $bin/reverse.pl $pro\.3UTR.up_minus.fa $pro\.3UTR.up_minus_rev.fa
cat $pro\.3UTR.up_plus.fa $pro\.3UTR.up_minus_rev.fa > $pro\.3UTR.up.fa
perl $bin/reverse.pl $pro\.3UTR.down_minus.fa $pro\.3UTR.down_minus_rev.fa
cat $pro\.3UTR.down_plus.fa $pro\.3UTR.down_minus_rev.fa > $pro\.3UTR.down.fa
perl /share/pub/xingsl/shilai/pipeline/tools/SignalSleuth2/SignalSleuth2.pl -seqfile $pro\.3UTR.up.fa -from 1 -to 100 -k 1 -gap 0 -sort T -cnt T
perl /share/pub/xingsl/shilai/pipeline/tools/SignalSleuth2/SignalSleuth2.pl -seqfile $pro\.3UTR.down.fa -from 1 -to 49 -k 1 -gap 0 -sort T -cnt T
sed -i 's/T/U/g' $pro\.3UTR.up.fa_1to100_k1_sort.cnt
sed -i 's/T/U/g' $pro\.3UTR.down.fa_1to49_k1_sort.cnt
/share/apps/R-4.2.1/bin/Rscript $bin/ATGC_plot.r $pro\.3UTR.up.fa_1to100_k1_sort.cnt $pro\.3UTR.down.fa_1to49_k1_sort.cnt $pro\.base.txt $pro\.base.png $pro\.base.pdf
rm $pro\.3UTR.up_plus.bed $pro\.3UTR.down_plus.bed $pro\.3UTR.up_minus.bed $pro\.3UTR.down_minus.bed 
rm $pro\.3UTR.up_plus.fa $pro\.3UTR.up_minus_rev.fa $pro\.3UTR.down_plus.fa $pro\.3UTR.down_minus.fa $pro\.3UTR.up_minus.fa $pro\.3UTR.down_minus_rev.fa 
#NUE
perl $bin/PA.up_down_strand.NUE.pl $pro\.intersect.3UTR.bed  $pro\.intersect.3UTR.NUE.bed
grep "+" $pro\.intersect.3UTR.NUE.bed  |cut -f1-3|sort -u > $pro\.intersect.3UTR.plus.NUE.bed
grep "-" $pro\.intersect.3UTR.NUE.bed  |cut -f1-3|sort -u > $pro\.intersect.3UTR.minus.NUE.bed
$bedtools getfasta -fi /MZAPA/Ref/B73_V4/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.intersect.3UTR.plus.NUE.bed > $pro\.intersect.3UTR.plus.NUE.fa
$bedtools getfasta -fi /MZAPA/Ref/B73_V4/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.intersect.3UTR.minus.NUE.bed >$pro\.intersect.3UTR.minus.NUE.fa
perl $bin/reverse.pl $pro\.intersect.3UTR.minus.NUE.fa $pro\.intersect.3UTR.minus_rev.NUE.fa
cat $pro\.intersect.3UTR.plus.NUE.fa $pro\.intersect.3UTR.minus_rev.NUE.fa >$pro\.intersect.3UTR.NUE.fa
rm $pro\.intersect.3UTR.plus.NUE.bed $pro\.intersect.3UTR.minus.NUE.bed
rm $pro\.intersect.3UTR.plus.NUE.fa $pro\.intersect.3UTR.minus_rev.NUE.fa $pro\.intersect.3UTR.minus.NUE.fa

