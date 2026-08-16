#python3 integrate_coloc_data.py --tissue Kern --qtl_type 3aeQTL --ld /APA/AFreview/coloc/S1960_vs_wanghaiyangNG/aeQTL_coloc/purity0.8/Kern/output_500kb_ld_cis_significant.txt --cis /APA/AFreview/coloc/S1960_vs_wanghaiyangNG/aeQTL_coloc/purity0.8/Kern/output_500kb.cis.tsv --anno /APA/AFreview/PAC_3seq_hq/region_dis_addStrand/SRP115041_ab22.intersect.3UTR.annoGene --output Kern.integrate_coloc.result.txt

for i in GRoot GShoot Kern L3Base L3Tip LMAD LMAN
do
	python3 integrate_coloc_data.py --tissue ${i} --qtl_type 3aeQTL --ld /APA/AFreview/coloc/S1960_vs_wanghaiyangNG/aeQTL_coloc/purity0.8/${i}/output_500kb_ld_cis_significant.txt --cis /APA/AFreview/coloc/S1960_vs_wanghaiyangNG/aeQTL_coloc/purity0.8/${i}/output_500kb.cis.tsv --anno /APA/AFreview/PAC_3seq_hq/region_dis_addStrand/SRP115041_ab22.intersect.3UTR.annoGene --output ${i}.integrate_coloc.result.txt
done
