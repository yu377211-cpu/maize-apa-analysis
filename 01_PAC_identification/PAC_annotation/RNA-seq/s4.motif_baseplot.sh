#!/bin/bash

bin=/APA/PAC_3seq_hq/motif/bin
for pro in ERP009123 ERP011069 ERP012784 Huang2018 Li2019 Lin2017 PRJNA171684_combine PRJNA189400 PRJNA277023 PRJNA437324 PRJNA482146 PRJNA505095 SRA055066 SRP006965 SRP017122 SRP018088 SRP026161 SRP037559 SRP044293 SRP069080
do
	dir=/MZAPA/ZY/post_analysis/PAC_CS_RNAseq_feature/$pro/filter01/DNA_motif/
	cd $dir
	#mv $pro\.base.png $pro\.base0.png
	#mv $pro\.base.pdf $pro\.base0.pdf
	/share/apps/R-4.2.1/bin/Rscript $bin/ATGC_plot.r $pro\.3UTR.up.fa_1to100_k1_sort.cnt $pro\.3UTR.down.fa_1to49_k1_sort.cnt $pro\.CS.base.txt $pro\.base.png $pro\.base.pdf
	#mv $pro\.CS.base.png $pro\.base.png
	#mv $pro\.CS.base.pdf $pro\.base.pdf
done
