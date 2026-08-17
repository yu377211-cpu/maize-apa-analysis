#!/bin/sh
#PBS -N map
#PBS -q workq
#PBS -l walltime=1000000:00:00
#PBS -l ncpus=4
#PBS -j oe
#PBS -m ae

sample=SRR5911572
seqdir=/data_download/MZAPA/SRP115041
outdir=/APA/PAC_3seq_map/map

mkdir $outdir/$sample\_2
cd $outdir/$sample\_2

#java -jar /share/apps/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 -threads 4 -summary summaryFile $seqdir/$sample/$sample\.fastq.gz $outdir/$sample/$sample\.clean.fastq ILLUMINACLIP:/share/apps/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:12 TRAILING:3 MINLEN:25
#perl /APA/PAC_3seq_map/PAC_identification_scripts/MAP_findTailAT.pl -in $outdir/$sample/$sample\.clean.fastq -poly 'A|T' -ml 25 -mp 8 -mg 5 -mm 2 -mr 2  -mtail 6 -debug F -odir $outdir/$sample/ -suf AT -mper 0.9 >$outdir/$sample/$sample\.findTailAT.o 2>$outdir/$sample/$sample\.findTailAT.e
#cat $outdir/$sample/$sample\.clean.AT.A.fq $outdir/$sample/$sample\.clean.AT.T.fq >$outdir/$sample/$sample\.clean.AT.fq
/share/apps/hisat2-2.1.0/hisat2 -p 4 -t --rna-strandness F -x /share2/pub/xingsl/xingsl/data_download/MZAPA/Ref/B73_V4/B73_V4 -U $outdir/$sample/$sample\.clean.AT.fq 2>$outdir/$sample\_2/$sample\.maplog -S $outdir/$sample\_2/$sample\.sam
#链特异性 如果使用dUTP: --rna-strandness R
/share/apps/samtools-1.9/bin/samtools view -F 20 $outdir/$sample\_2/$sample\.sam |grep "NH:i:1" >$outdir/$sample\_2/$sample\.F.sam
/share/apps/samtools-1.9/bin/samtools view -f 16 $outdir/$sample\_2/$sample\.sam |grep "NH:i:1" >$outdir/$sample\_2/$sample\.R.sam
cat $outdir/$sample\_2/$sample\.F.sam $outdir/$sample\_2/$sample\.R.sam >$outdir/$sample\_2/$sample\.umap.sam
#rm $outdir/$sample/$sample\.F.sam  $outdir/$sample/$sample\.R.sam
#rm $outdir/$sample/*.raw.fq
#ln -sf $outdir/$sample/$sample\.umap.sam $outdir/sam_file/

