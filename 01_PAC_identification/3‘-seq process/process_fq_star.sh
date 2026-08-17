#!/bin/sh
#PBS -N map
#PBS -q workq
#PBS -l walltime=1000000:00:00
#PBS -l ncpus=6
#PBS -j oe
#PBS -m ae

sample=${sample}
seqdir=/data_download/MZAPA/SRP115041
outdir=/APA/PAC_3seq_map/map

mkdir $outdir/$sample
cd $outdir/$sample

java -jar /share/apps/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 -threads 4 -summary summaryFile $seqdir/$sample/$sample\.fastq.gz $outdir/$sample/$sample\.clean.fastq ILLUMINACLIP:/share/apps/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:12 TRAILING:3 MINLEN:25
perl /APA/PAC_3seq_map/PAC_identification_scripts/MAP_findTailAT.pl -in $outdir/$sample/$sample\.clean.fastq -poly 'A|T' -ml 25 -mp 8 -mg 5 -mm 2 -mr 2  -mtail 6 -debug F -odir $outdir/$sample/ -suf AT -mper 0.9 >$outdir/$sample/$sample\.findTailAT.o 2>$outdir/$sample/$sample\.findTailAT.e
cat $outdir/$sample/$sample\.clean.AT.A.fq $outdir/$sample/$sample\.clean.AT.T.fq >$outdir/$sample/$sample\.clean.AT.fq
/share/apps/STAR-2.6.0a/bin/Linux_x86_64/STAR --runThreadN 12 \
     --readFilesIn $outdir/$sample/$sample\.clean.AT.fq \
     --genomeDir /share2/pub/xingsl/xingsl/data_download/MZAPA/Ref/B73_V4/ \
     --outFileNamePrefix   $outdir/$sample/$sample \
     --outMultimapperOrder Random \
     --outFilterMultimapNmax 1
rm  $outdir/$sample/$sample\.clean.AT.A.fq $outdir/$sample/$sample\.clean.AT.T.fq
