bin=/home/zhangyu/workdir/APA/post_analysis/3seq/SRP115041/bin
grep Chr /home/zhangyu/workdir/APA/post_analysis/3seq/SRP115041/PAC_umap_all_A/all.PAC.PATcount >umap.Chr.PAC.PATcount
cat /home/zhangyu/workdir/APA/post_analysis/3seq/SRP115041/PAC_umap_all_A/all.PAC.header umap.Chr.PAC.PATcount|awk '{print "PAC"NR-1"\t"$0}' >umap.Chr.PAC.PATcount.txt
awk '{print "PAC"NR"\t"$3-$2+1}' umap.Chr.PAC.PATcount >umap.Chr.PAC.length
awk -F"\t" '{for (i=10;i<=NF;i++)printf("%s\t", $i);print ""}' umap.Chr.PAC.PATcount.txt|awk '{print "PAC"NR-1"\t"$0}' >umap.Chr.PAC.PATcount.matrix
sed -i s/PAC0/PACid/g umap.Chr.PAC.PATcount.matrix
sed -i 's/[\t ]\+$//' umap.Chr.PAC.PATcount.matrix
perl $bin/cal_TPM.pl umap.Chr.PAC.length umap.Chr.PAC.PATcount.matrix umap.Chr.TPM
sed -i 's/[\t ]\+$//' umap.Chr.TPM
#计算PAC的平均TPM
awk 'BEGIN{sum = 0} {for(i = 2; i <= NF; i++) {sum += $i} {print $1"\t"sum/(NF-1); sum = 0}}' umap.Chr.TPM >umap.Chr.mean.TPM
#PAC样本表达比
awk 'BEGIN{sum = 0} {for(i = 2; i <= NF; i++) {if($i >=2) sum += 1} {print $1"\t"sum"/1960"; sum = 0}}' umap.Chr.PAC.PATcount.matrix >umap.Chr.SampleExpressRatio
#高质量PAC
awk 'BEGIN{sum = 0} {for(i = 2; i <= NF; i++) {if($i >=1) sum += 1} {if(sum >=2 )print $0; sum = 0}}' umap.Chr.TPM >umap.Chr.hq.TPM
perl $bin/tiqu.pl umap.Chr.hq.TPM umap.Chr.PAC.PATcount.txt umap.Chr.PAC.hq.pos
#方向一致PAC
umap.Chr.PAC.hq.pos ==》umap.Chr.hq.bed
perl tiqu2.pl umap.Chr.hq.bed /home/zhangyu/workdir/APA/post_analysis/3seq/SRP115041/region_dis_A/transcript_new/SRP115041.intersect.alltype.ori_cons umap.Chr.hq.ori_cons.bed
perl tiqu3.pl umap.Chr.hq.ori_cons.bed umap.Chr.hq.TPM umap.Chr.hq.ori_cons.TPM
perl add_info.pl /home/zhangyu/workdir/APA/post_analysis/3seq/SRP115041/region_dis_A/transcript_new/SRP115041.intersect.alltype.uniq umap.Chr.hq.bed umap.Chr.hq.bed.anno

cut -f1,5 umap.Chr.hq.bed >umap.Chr.hq.pos
python3 /home/zhangyu/workdir/APA/post_analysis/GWAS/site_annotation.py umap.Chr.hq.pos /home/zhangyu/workdir/APA/post_analysis/GWAS/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gene.bed umap.Chr.hq.pos.anno
cut -f1,5 umap.Chr.hq.ori_cons.bed >umap.Chr.hq.ori_cons.pos
python3 /home/zhangyu/workdir/APA/post_analysis/GWAS/site_annotation.py umap.Chr.hq.ori_cons.pos /home/zhangyu/workdir/APA/post_analysis/GWAS/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gene.bed umap.Chr.hq.ori_cons.pos.anno

sh pac_of_samples.sh umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.txt
perl sample_base_pac.pl umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.txt.sampleSta SRP115041.ReadsQuality.xls SRP115041.sample_base.txt
/share/apps/R-4.2.1/bin/Rscript corelation_canshu.r

awk 'BEGIN{sum = 0} {for(i = 2; i <= NF; i++) {if($i > 0) sum += 1} {print $1"\t"sum"\t1960"; sum = 0}}' umap.Chr.PAC.PATcount.matrix >umap.Chr.SampleExpressRatio_ab0
awk -F"\t" '{ sum += $2; n++ } END { print "Average:", sum / n }' umap.Chr.SampleExpressRatio_ab0
Average: 21.5317
awk -F"\t" '{if($2 >= 22){print $1}}' umap.Chr.SampleExpressRatio_ab0 >umap.Chr.SampleExpressRatio_ab22.id
perl tiqu_PAC.pl umap.Chr.SampleExpressRatio_ab22.id umap.Chr.PAC.PATcount.txt umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.txt
wc -l umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.txt
168383
awk 'BEGIN{sum = 0} {for(i = 10; i <= NF; i++) {if($i >=1) sum += 1} {if(sum >=2 )print $0; sum = 0}}' umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.txt|wc -l
168383
perl tiqu_PAC.pl umap.Chr.SampleExpressRatio_ab22.id umap.Chr.hq.TPM umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.TPM
cut -f1,10- umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.txt >umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.matrix
perl cal_TPM.pl umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.len umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.matrix umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.TPM2
sed -i 's/[\t ]\+$//' umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.TPM2
