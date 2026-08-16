#!/bin/bash

#GRoot GShoot Kern L3Base L3Tip LMAD LMAN
for name in GShoot
#for name in GRoot Kern L3Base L3Tip LMAD LMAN
do
	VCF="/APA/3aQTL_link/3seq_PDUI_TPM/${name}_addM/phenotype/S296_v4.${name}.vcf.gz"
	dir="/APA/AFreview/coloc/S1960_vs_wanghaiyangNG/aeQTL_coloc"
	OUTPUT=${dir}/${name}
	bcftools="/share/apps/bcftools-1.9/bin/bcftools"

	cd $OUTPUT

	cut -f1,2,9 $OUTPUT/output_500kb.cis.tsv|awk -F"\t" '{print $1"_"$2"\t"$1"_"$3}' |grep -v "chr.x"|sort -u > ${OUTPUT}/output_500kb.cis.pair
	cat output_500kb.cis.pair | tr '\t' '\n' | sort | uniq >all_cis_target_snps.txt
	$bcftools view -i 'ID=@all_cis_target_snps.txt' "$VCF" -Oz -o "${OUTPUT}.cis.vcf.gz"
	$bcftools index "${OUTPUT}.cis.vcf.gz"
	python3 /APA/3aQTL_link/coloc/bin/ld_calculator.py  --vcf ${OUTPUT}.cis.vcf.gz --pairs ${OUTPUT}/output_500kb.cis.pair --output ${OUTPUT}/output_500kb_ld_cis_results.tsv --plink /share/pub/xingsl/shilai/pipeline/tools/plink_1.9_20250819/plink --processes 8 --tmpdir ${OUTPUT}/ld_logs_cis
	awk -F"\t" '{if($3 >0.8){print}}' ${OUTPUT}/output_500kb_ld_cis_results.tsv|grep -v NA >${OUTPUT}/output_500kb_ld_cis_significant.txt
done

