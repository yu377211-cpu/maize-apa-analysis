#!/bin/bash

#GRoot GShoot Kern L3Base L3Tip LMAD LMAN
for name in GShoot

#for name in GRoot Kern L3Base L3Tip LMAD LMAN
do
	myQTL="/APA/AFreview/MeQTL/${name}/FineMapping/susieR_res.all_genes.txt"
	paperdata="/APA/3aQTL_link/coloc/15.S1960_vs_wanghaiyangNG/NG_site.v4.format.txt"
	mkdir /APA/AFreview/coloc/S1960_vs_wanghaiyangNG/aeQTL_coloc/$name
	dir="/APA/AFreview/coloc/S1960_vs_wanghaiyangNG/aeQTL_coloc/$name"

	#awk -F'\t' 'BEGIN {OFS="\t"; print "chr", "pos", "locus_id", "variant_id", "pip", "cs", "cs_size", "cs_purity"} NR==1 {next} {split($2,a, "_"); print a[1], a[2], $1, $2, $3, $4, $5, $6}' $myQTL > ${name}.susieR_res.all_genes.txt
	awk -F"\t" '{if($6 >=0.8){print}}' $myQTL |grep -v NA|awk -F'\t' 'BEGIN {OFS="\t"; print "chr", "pos", "locus_id", "variant_id", "pip", "cs", "cs_size", "cs_purity"} NR==1 {next} {split($2,a, "_"); print a[1], a[2], $1, $2, $3, $4, $5, $6}' >${name}.susieR_res.filter_genes.txt
	Rscript /APA/3aQTL_link/coloc/bin/find_500kb.r ${name}.susieR_res.filter_genes.txt $paperdata $dir/output_500kb.cis.tsv 500000
done

#cut -f5 $dir/output_500kb.cis.tsv |sort -g|awk '{a[NR]=$1} END{if(NR%2==1) print "中位数: " a[(NR+1)/2]; else print "中位数: " (a[NR/2]+a[NR/2+1])/2}'
#中位数: 2.15e-18
#cut -f8 $dir/output_500kb.cis.tsv |sort -g|awk '{a[NR]=$1} END{if(NR%2==1) print "中位数: " a[(NR+1)/2]; else print "中位数: " (a[NR/2]+a[NR/2+1])/2}'
#中位数: 0.02
#Rscript /APA/3aQTL_link/coloc/bin/genome_coloc_canshu.r $dir/output_500kb.cis.tsv /APA/3aQTL_link/coloc/bin/chrNameLength_sub.txt $dir/Seedling_vs_${name}_cis_genome_wide_plot
#中位数: 1.5148284905567498e-18
#中位数: 0.024531129

