pro=PRJNA277023
bin=/home/zhangyu/workdir/APA/post_analysis/compare/bin

mkdir motif
cd motif
##CS
grep "three_prime_UTR" ../$pro\.intersect.alltype|awk -F"\t" '{OFS="\t";print $1,$2,$3,$4,".",$12,$10}'|awk -F"_" '{print $1}'|sort -u >$pro\.PAC.CS.3UTR.gene_strand.bed 
perl /home/zhangyu/workdir/APA/post_analysis/compare/bin/PA.up_down_strand.pl $pro\.PAC.CS.3UTR.gene_strand.bed $pro\.CS_3UTR.ud.bed
grep "+" $pro\.CS_3UTR.ud.bed|grep upstream |cut -f1-3|sort -u >$pro\.CS_3UTR.up_plus.bed
grep "+" $pro\.CS_3UTR.ud.bed|grep downstream|cut -f1-3|sort -u >$pro\.CS_3UTR.down_plus.bed
grep "-" $pro\.CS_3UTR.ud.bed|grep upstream |cut -f1-3|sort -u >$pro\.CS_3UTR.up_minus.bed
grep "-" $pro\.CS_3UTR.ud.bed|grep downstream|cut -f1-3|sort -u >$pro\.CS_3UTR.down_minus.bed

/public/software/bedtools2/bin/bedtools getfasta -fi /home/zhangyu/workdir/APA/maize_ref/fromXSL/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.CS_3UTR.up_plus.bed >$pro\.CS_3UTR.up_plus.fa
/public/software/bedtools2/bin/bedtools getfasta -fi /home/zhangyu/workdir/APA/maize_ref/fromXSL/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.CS_3UTR.down_plus.bed >$pro\.CS_3UTR.down_plus.fa
/public/software/bedtools2/bin/bedtools getfasta -fi /home/zhangyu/workdir/APA/maize_ref/fromXSL/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.CS_3UTR.up_minus.bed >$pro\.CS_3UTR.up_minus.fa
/public/software/bedtools2/bin/bedtools getfasta -fi /home/zhangyu/workdir/APA/maize_ref/fromXSL/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.CS_3UTR.down_minus.bed >$pro\.CS_3UTR.down_minus.fa
perl /home/zhangyu/workdir/APA/post_analysis/compare/bin/reverse.pl $pro\.CS_3UTR.up_minus.fa $pro\.CS_3UTR.up_minus_rev.fa
cat $pro\.CS_3UTR.up_plus.fa $pro\.CS_3UTR.up_minus_rev.fa > $pro\.CS_3UTR.up.fa
perl /home/zhangyu/workdir/APA/post_analysis/compare/bin/reverse.pl $pro\.CS_3UTR.down_minus.fa $pro\.CS_3UTR.down_minus_rev.fa
cat $pro\.CS_3UTR.down_plus.fa $pro\.CS_3UTR.down_minus_rev.fa > $pro\.CS_3UTR.down.fa
perl /home/zhangyu/workdir/software/SignalSleuth2/SignalSleuth2.pl -seqfile $pro\.CS_3UTR.up.fa -from 1 -to 100 -k 1 -gap 0 -sort T -cnt T 
perl /home/zhangyu/workdir/software/SignalSleuth2/SignalSleuth2.pl -seqfile $pro\.CS_3UTR.down.fa -from 1 -to 49 -k 1 -gap 0 -sort T -cnt T
sed -i 's/T/U/g' $pro\.CS_3UTR.up.fa_1to100_k1_sort.cnt
sed -i 's/T/U/g' $pro\.CS_3UTR.down.fa_1to49_k1_sort.cnt
Rscript $bin/ATGC_plot.r $pro\.CS_3UTR.up.fa_1to100_k1_sort.cnt $pro\.CS_3UTR.down.fa_1to49_k1_sort.cnt $pro\.CS.base.txt $pro\.CS.base.png $pro\.CS.base.pdf

##DP2
awk -F"\t" '{if($7 != "."){OFS="\t";print $1,$2,$3,$4,$5,$12,$10}}' ../$pro\_DP2.intersect.3UTR|awk -F"_" '{print $1"_"$2}'|sort -u >$pro\.PAC.DP2.3UTR.gene_strand.bed
perl /home/zhangyu/workdir/APA/post_analysis/compare/bin/PA.up_down_strand.pl $pro\.PAC.DP2.3UTR.gene_strand.bed $pro\.DP2_PA3UTR.ud.bed
grep "+" $pro\.DP2_PA3UTR.ud.bed|grep upstream |cut -f1-3|sort -u > $pro\.DP2_PA3UTR.up_plus.bed
grep "+" $pro\.DP2_PA3UTR.ud.bed|grep downstream|cut -f1-3|sort -u > $pro\.DP2_PA3UTR.down_plus.bed
grep "-" $pro\.DP2_PA3UTR.ud.bed|grep upstream |cut -f1-3|sort -u > $pro\.DP2_PA3UTR.up_minus.bed
grep "-" $pro\.DP2_PA3UTR.ud.bed|grep downstream |cut -f1-3|sort -u > $pro\.DP2_PA3UTR.down_minus.bed

/public/software/bedtools2/bin/bedtools getfasta -fi /home/zhangyu/workdir/APA/maize_ref/fromXSL/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.DP2_PA3UTR.up_plus.bed >$pro\.DP2_PA3UTR.up_plus.fa
/public/software/bedtools2/bin/bedtools getfasta -fi /home/zhangyu/workdir/APA/maize_ref/fromXSL/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.DP2_PA3UTR.down_plus.bed >$pro\.DP2_PA3UTR.down_plus.fa
/public/software/bedtools2/bin/bedtools getfasta -fi /home/zhangyu/workdir/APA/maize_ref/fromXSL/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.DP2_PA3UTR.up_minus.bed >$pro\.DP2_PA3UTR.up_minus.fa
/public/software/bedtools2/bin/bedtools getfasta -fi /home/zhangyu/workdir/APA/maize_ref/fromXSL/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.DP2_PA3UTR.down_minus.bed >$pro\.DP2_PA3UTR.down_minus.fa

perl /home/zhangyu/workdir/APA/post_analysis/compare/bin/reverse.pl $pro\.DP2_PA3UTR.up_minus.fa $pro\.DP2_PA3UTR.up_minus_rev.fa
cat $pro\.DP2_PA3UTR.up_plus.fa $pro\.DP2_PA3UTR.up_minus_rev.fa >$pro\.DP2_PA3UTR.up.fa
perl /home/zhangyu/workdir/APA/post_analysis/compare/bin/reverse.pl $pro\.DP2_PA3UTR.down_minus.fa $pro\.DP2_PA3UTR.down_minus_rev.fa
cat $pro\.DP2_PA3UTR.down_plus.fa $pro\.DP2_PA3UTR.down_minus_rev.fa >$pro\.DP2_PA3UTR.down.fa
perl /home/zhangyu/workdir/software/SignalSleuth2/SignalSleuth2.pl -seqfile $pro\.DP2_PA3UTR.up.fa -from 1 -to 100 -k 1 -gap 0 -sort T -cnt T
perl /home/zhangyu/workdir/software/SignalSleuth2/SignalSleuth2.pl -seqfile $pro\.DP2_PA3UTR.down.fa -from 1 -to 49 -k 1 -gap 0 -sort T -cnt T
sed -i 's/T/U/g' $pro\.DP2_PA3UTR.up.fa_1to100_k1_sort.cnt
sed -i 's/T/U/g' $pro\.DP2_PA3UTR.down.fa_1to49_k1_sort.cnt
Rscript $bin/ATGC_plot.r $pro\.DP2_PA3UTR.up.fa_1to100_k1_sort.cnt $pro\.DP2_PA3UTR.down.fa_1to49_k1_sort.cnt $pro\.DP2.base.txt $pro\.DP2.base.png $pro\.DP2.base.pdf

#MEME
#DP2
<<BLOCK
perl /home/zhangyu/workdir/APA/post_analysis/compare/bin/PA.up_down_strand.2.pl $pro\.PAC.DP2.3UTR.gene_strand.bed $pro\.PAC.DP2.3UTR.motif.bed
grep "+" $pro\.PAC.DP2.3UTR.motif.bed|cut -f1-3|sort -u > $pro\.PAC.DP2.3UTR.motif_plus.bed
grep "-" $pro\.PAC.DP2.3UTR.motif.bed|cut -f1-3|sort -u > $pro\.PAC.DP2.3UTR.motif_minus.bed
/public/software/bedtools2/bin/bedtools getfasta -fi /home/zhangyu/workdir/APA/maize_ref/fromXSL/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.PAC.DP2.3UTR.motif_plus.bed >$pro\.PAC.DP2.3UTR.motif_plus.fa
/public/software/bedtools2/bin/bedtools getfasta -fi /home/zhangyu/workdir/APA/maize_ref/fromXSL/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.PAC.DP2.3UTR.motif_minus.bed >$pro\.PAC.DP2.3UTR.motif_minus.fa
perl /home/zhangyu/workdir/APA/post_analysis/compare/bin/reverse.pl $pro\.PAC.DP2.3UTR.motif_minus.fa $pro\.PAC.DP2.3UTR.motif_minus_rev.fa
cat $pro\.PAC.DP2.3UTR.motif_plus.fa $pro\.PAC.DP2.3UTR.motif_minus_rev.fa >$pro\.PAC.DP2.3UTR.motif.fa
meme $pro\.PAC.DP2.3UTR.motif.fa -dna -oc DP2 -nostatus -time 14400 -mod zoops -nmotifs 10 -minw 5 -maxw 25 -objfun classic -revcomp -markov_order 0 -evt 0.05 -o DP2.meme
#CS meme
perl /home/zhangyu/workdir/APA/post_analysis/compare/bin/PA.up_down_strand.2.pl $pro\.PAC.CS.3UTR.gene_strand.bed $pro\.PAC.CS.3UTR.motif.bed
grep "+" $pro\.PAC.CS.3UTR.motif.bed |cut -f1-3|sort -u >$pro\.PAC.CS.3UTR.motif_plus.bed
grep "-" $pro\.PAC.CS.3UTR.motif.bed |cut -f1-3|sort -u >$pro\.PAC.CS.3UTR.motif_minus.bed
/public/software/bedtools2/bin/bedtools getfasta -fi /home/zhangyu/workdir/APA/maize_ref/fromXSL/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.PAC.CS.3UTR.motif_plus.bed >$pro\.PAC.CS.3UTR.motif_plus.fa
/public/software/bedtools2/bin/bedtools getfasta -fi /home/zhangyu/workdir/APA/maize_ref/fromXSL/Zm-B73-REFERENCE-GRAMENE-4.0.fa -bed $pro\.PAC.CS.3UTR.motif_minus.bed >$pro\.PAC.CS.3UTR.motif_minus.fa
perl /home/zhangyu/workdir/APA/post_analysis/compare/bin/reverse.pl $pro\.PAC.CS.3UTR.motif_minus.fa $pro\.PAC.CS.3UTR.motif_minus_rev.fa
cat $pro\.PAC.CS.3UTR.motif_plus.fa $pro\.PAC.CS.3UTR.motif_minus_rev.fa >$pro\.PAC.CS.3UTR.motif.fa
meme $pro\.PAC.CS.3UTR.motif.fa -dna -oc CS -nostatus -time 14400 -mod zoops -nmotifs 10 -minw 5 -maxw 25 -objfun classic -revcomp -markov_order 0 -evt 0.05 -o CS.meme
BLOCK
meme $pro\.DP2_PA3UTR.up.fa -dna -oc DP2 -nostatus -time 14400 -mod zoops -nmotifs 10 -minw 5 -maxw 25 -objfun classic -revcomp -markov_order 0 -evt 0.05 -o DP2.meme
meme $pro\.CS_3UTR.up.fa -dna -oc CS -nostatus -time 14400 -mod zoops -nmotifs 10 -minw 5 -maxw 25 -objfun classic -revcomp -markov_order 0 -evt 0.05 -o CS.meme
