pro=$1
dir=/MZAPA/ZY/post_analysis/PAC_CS_RNAseq_feature
cd $dir/$pro/filter01

awk -F"\t" '{if($7 != "."){OFS="\t";print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,"ncRNA."$11,$12}}' $pro\.intersect.ncRNA > $pro\.intersect.ncRNA.clean
cat $pro\.intersect.3UTR $pro\.intersect.5UTR $pro\.intersect.CDS $pro\.intersect.extend $pro\.intersect.intergenic $pro\.intersect.intron $pro\.intersect.ncRNA.clean |awk -F"\t" '{if($7 != "."){print}}' > $pro\.intersect.allType
cut -f5,11 $pro\.intersect.allType|sort -u >$pro\.intersect.allType.txt
sed -i s/\.lincRNA//g $pro\.intersect.allType.txt
sed -i s/\.miRNA//g $pro\.intersect.allType.txt
sed -i s/\.protein_coding//g $pro\.intersect.allType.txt
sed -i s/\.tRNA//g $pro\.intersect.allType.txt
sed -i s/\.transcript//g $pro\.intersect.allType.txt
sed -i s/five_prime_UTR/five_UTR/g $pro\.intersect.allType.txt
sed -i s/three_prime_UTR/three_UTR/g $pro\.intersect.allType.txt
sort -u $pro\.intersect.allType.txt >$pro\.intersect.allType.uniq.txt
rm $pro\.intersect.allType.txt

awk -F"\t" '{print $5"\t"$3-$2+1}' $dir/$pro/filter01/$pro\.hq.PAC > $dir/$pro/filter01/$pro\.hq.PAC.len 
sed -i '1d' $dir/$pro/filter01/$pro\.hq.PAC.len
sed -i '1i PACid\tlen' $dir/$pro/filter01/$pro\.hq.PAC.len

grep -v "intergenic"  $pro\.intersect.allType|cut -f5,10|awk -F"_" '{print $1}' |sort -u|cut -f2|sort |uniq -c|sort -k1n >$pro\.hq.gene.sta
echo "$pro"
awk '{if($1 == "1"){print}}' $pro\.hq.gene.sta|wc -l
awk '{if($1 == "2"){print}}' $pro\.hq.gene.sta|wc -l
awk '{if($1 > 2){print}}' $pro\.hq.gene.sta|wc -l

/share/apps/R-4.2.1/bin/Rscript /APA/PAC_3seq_hq/sta_plot/box_single_TPM.r $pro\.intersect.allType.uniq.txt $pro\.hq.PAC.TPM $pro\.PAC_TPM.box
/share/apps/R-4.2.1/bin/Rscript /APA/PAC_3seq_hq/sta_plot/box_single_len.r $pro\.intersect.allType.uniq.txt $pro\.hq.PAC.len $pro\.PAC_len.box
