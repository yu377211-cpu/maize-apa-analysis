pro=${pro}
indir=/MZAPA/ZY/CS_RNAseq/CS_all.DataSet
outdir=/MZAPA/ZY/post_analysis/PAC_CS_RNAseq_feature/$pro
data_dir=/MZAPA/Ref/B73_V4
bedtools=/share/apps/bedtools2/bin/bedtools
bin=/MZAPA/ZY/post_analysis/PAC_CS_RNAseq_feature/bin
mkdir $outdir
mkdir $outdir/filter01
cd $outdir/filter01
##filter 
grep Chr $indir/$pro\.all.PAC.PATcount.txt|awk -F"\t" '{if($9 > 2){OFS="\t";print }}' >$pro\.hq1.PATcount
#PAC1120 Chr1    47059   47059   +       1       45      47059   45
awk '{print $1"\t"$4-$3+1}' $pro\.hq1.PATcount > $pro\.hq1.PAC.length
head -n1 $indir/$pro\.all.PAC.PATcount.txt >title
cat title $pro\.hq1.PATcount  >$pro\.hq1.PATcount.txt
cut -f1,10- $pro\.hq1.PATcount.txt > $pro\.hq1.PATcount.matrix
sed -i s/PAC0/PACid/g $pro\.hq1.PATcount.matrix
perl $bin/cal_TPM.pl $pro\.hq1.PAC.length $pro\.hq1.PATcount.matrix $pro\.hq1.PAC.TPM
sed -i 's/[\t ]\+$//' $pro\.hq1.PAC.TPM
cd $outdir/filter01
awk 'BEGIN{sum = 0} {for(i = 2; i <= NF; i++) {if($i >=1) sum += 1} {if(sum >=2 )print $0; sum = 0}}' $pro\.hq1.PAC.TPM > $pro\.hq.PAC.TPM
perl $bin/tiqu_PAC.pl $pro\.hq.PAC.TPM  $indir/$pro\.all.PAC.PATcount.txt $pro\.hq.PAC
sort -V $pro\.hq.PAC >$pro\.hq.PAC.s
mv $pro\.hq.PAC.s $pro\.hq.PAC
rm $pro\.hq1.*
rm title

##3’UTR>extend>CDS>intron>5’UTR>ncRNA>intergenic
$bedtools intersect -a $pro\.hq.PAC  -b $data_dir/three_prime_UTR.bed -wa -wb -loj >$pro\.intersect.3UTR
awk -F"\t" '{if($7 == "."){print}}' $pro\.intersect.3UTR|cut -f1-6 >round1.pos
$bedtools intersect -a round1.pos -b $data_dir/extend.bed -wa -wb -loj >$pro\.intersect.extend
awk -F"\t" '{if($7 == "."){print}}' $pro\.intersect.extend |cut -f1-6 >round2.pos
$bedtools intersect -a round2.pos -b $data_dir/CDS.bed -wa -wb -loj >$pro\.intersect.CDS
awk -F"\t" '{if($7 == "."){print}}' $pro\.intersect.CDS |cut -f1-6 >round3.pos
$bedtools intersect -a round3.pos -b $data_dir/intron.bed -wa -wb -loj >$pro\.intersect.intron
awk -F"\t" '{if($7 == "."){print}}' $pro\.intersect.intron |cut -f1-6 >round4.pos
$bedtools intersect -a round4.pos -b $data_dir/five_prime_UTR.bed -wa -wb -loj >$pro\.intersect.5UTR
awk -F"\t" '{if($7 == "."){print}}' $pro\.intersect.5UTR|cut -f1-6 >round5.pos
$bedtools intersect -a round5.pos -b $data_dir/ncRNA.bed -wa -wb -loj >$pro\.intersect.ncRNA
awk -F"\t" '{if($7 == "."){print}}' $pro\.intersect.ncRNA|cut -f1-6 >round6.pos
$bedtools intersect -a round6.pos -b $data_dir/intergenic_addinfo.bed -wa -wb -loj >$pro\.intersect.intergenic

#统计各类PAC数目
echo $pro
for i in { 3UTR extend CDS intron 5UTR ncRNA intergenic }
do
	awk  -F"\t" '{if($7 != "."){print}}' $pro\.intersect.$i|cut -f1-4|sort -u|wc -l
done

for i in 3UTR extend CDS intron 5UTR ncRNA intergenic
do
	awk  -F"\t" '{if($7 != "."){print}}' $pro\.intersect.$i|cut -f1-6|sort -u >$pro\.intersect.$i\.clean
done
