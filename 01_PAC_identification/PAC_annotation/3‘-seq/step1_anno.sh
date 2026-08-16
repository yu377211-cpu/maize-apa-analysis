#!/bin/sh

bedtools=/share/apps/bedtools2/bin/bedtools
dir=/APA/PAC_3seq_hq/region_dis

info= "/APA/PAC_3seq_hq/umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.txt"
awk -F"\t" '{OFS="\t"; print $2,$3,$4,$8,$1,$5}' /APA/PAC_3seq_hq/umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.txt >umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.txt.pos
sort -V umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.txt.pos >umap.Chr.PAC.SampleExpressRatio_ab22.sort.pos
rm umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.txt.pos
ln -s /APA/PAC_3seq_hq/region_dis/umap.Chr.PAC.SampleExpressRatio_ab22.sort.pos ./
##3’UTR>extend>CDS>intron>5’UTR>ncRNA>intergenic
$bedtools intersect -s -a umap.Chr.PAC.SampleExpressRatio_ab22.sort.pos  -b $dir/three_prime_UTR.bed -wa -wb -loj >SRP115041_ab22.intersect.3UTR
awk -F"\t" '{if($7 == "."){print}}' SRP115041_ab22.intersect.3UTR|cut -f1-6 >round1.pos
/share/apps/bedtools2.25.0/bin/bedtools intersect -s -a round1.pos -b $dir/extend.bed -wa -wb -loj >SRP115041_ab22.intersect.extend
awk -F"\t" '{if($7 == "."){print}}' SRP115041_ab22.intersect.extend |cut -f1-6 >round2.pos
/share/apps/bedtools2.25.0/bin/bedtools intersect -s -a round2.pos -b $dir/CDS.bed -wa -wb -loj >SRP115041_ab22.intersect.CDS
awk -F"\t" '{if($7 == "."){print}}' SRP115041_ab22.intersect.CDS |cut -f1-6 >round3.pos
/share/apps/bedtools2.25.0/bin/bedtools intersect -s -a round3.pos -b $dir/intron.bed -wa -wb -loj >SRP115041_ab22.intersect.intron
awk -F"\t" '{if($7 == "."){print}}' SRP115041_ab22.intersect.intron |cut -f1-6 >round4.pos
/share/apps/bedtools2.25.0/bin/bedtools intersect -s -a round4.pos -b $dir/five_prime_UTR.bed -wa -wb -loj >SRP115041_ab22.intersect.5UTR
awk -F"\t" '{if($7 == "."){print}}' SRP115041_ab22.intersect.5UTR|cut -f1-6 >round5.pos
/share/apps/bedtools2.25.0/bin/bedtools intersect -s -a round5.pos -b $dir/ncRNA.bed -wa -wb -loj >SRP115041_ab22.intersect.ncRNA
awk -F"\t" '{if($7 == "."){print}}' SRP115041_ab22.intersect.ncRNA|cut -f1-6 >round6.pos
/share/apps/bedtools2.25.0/bin/bedtools intersect -a round6.pos -b $dir/intergenic_addinfo.bed -wa -wb -loj >SRP115041_ab22.intersect.intergenic

for i in 3UTR extend CDS intron 5UTR ncRNA intergenic ;
do
	awk -F"\t" '{if($7 != "."){print}}' SRP115041_ab22.intersect.$i|cut -f1-4|sort -u|wc -l
done


