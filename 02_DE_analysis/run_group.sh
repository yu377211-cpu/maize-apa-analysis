#!/bin/bash


dir=/APA/AFreview/PAC_3seq_diffana/PAC_readcount/same_tissue/diff
bin=/APA/PAC_3seq_diffana/bin
group=" GRoot GShoot Kern L3Base L3Tip LMAD LMAN "
lst=" mixedvsnss mixedvspopcorn mixedvsss mixedvssweet mixedvsts nssvspopcorn nssvsss nssvssweet nssvsts popcornvsss popcornvssweet popcornvsts ssvssweet ssvsts sweetvsts lnvsmixed lnvsnss lnvspopcorn lnvsss lnvssweet lnvsts "

cd $dir || exit

for subset in ${group}
do
	cd $dir || exit
#	grep -f ${dir}/../${subset}.lst /APA/PAC_3seq_diffana/S1960.group.txt >${subset}.group
#	ln -s /APA/PAC_3seq_diffana/PAC_tpm_ab22/same_tissue/diff/${subset}.group ./
	mkdir -p $subset
	cd $dir/$subset

#	for i in ${lst}
#	do
#		python3 $bin/get_group_condition.py -i ${dir}/${subset}.group -c $i -o ${i}.condition
#	done

	export R_LIBS="/share/apps/R-3.6.3/lib64/R/library/":$R_LIBS
	export R_LIBS="/APA/software/miniforge3/envs/R3.6.3/lib/R/library/":$R_LIBS
	for i in ${lst}
	do
		/APA/software/miniforge3/envs/R3.6.3/bin/Rscript $bin/DESeq2_covariate.R $dir/../${subset}.hq.3UTR.tag ${i}.condition ${subset}.${i}.heatmap.png "~ group"
	done
done
