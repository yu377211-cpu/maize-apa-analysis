#!/bin/sh
#PBS -N aeQTL_GRoot
#PBS -q workq
#PBS -l walltime=1000000:00:00
#PBS -l ncpus=12
#PBS -j oe
#PBS -m ae

group=GRoot
out_dir=/APA/AFreview/MeQTL
cpm=/APA/AFreview/PAC_3seq_hq/region_dis_addStrand/umap.Chr.ab22.matrix_TMM_CPM.3UTR.txt
bin=/APA/AFreview/3aQTL/bin
<<BLOCK
mkdir -p $out_dir/$group
mkdir -p $out_dir/$group/phenotype
cd  $out_dir/$group/phenotype
python3 /APA/PAC_3seq_diffana/bin/find_row2.py -k /APA/MeQTL/union_1960_AGP/$group\.final.comb.lst -f $cpm -o $group\.paired.cpm
perl $bin/PAC.changeID.pl $group\.paired.cpm /APA/MeQTL/union_1960_AGP/$group\.final.comb $group\.paired.cpm
sed -i 's/\t$//g' $group\.paired.cpm.changeID.txt
head -n1 /APA/MeQTL/genotype/$group\.Genotype.eQTL.2.txt|tr '\t' '\n' > $out_dir/$group/phenotype/$group\.lst
ln -s /APA/3aQTL_link/3seq_PDUI_TPM/${group}_addM/phenotype/S296_v4.${group}.vcf.gz ./
ln -s /APA/3aQTL_link/3seq_PDUI_TPM/${group}_addM/phenotype/S296_v4.${group}.vcf.gz.tbi ./
echo "$out_dir/$group/phenotype/S296_v4.${group}.vcf.gz" >$out_dir/$group/vcf_list.txt
BLOCK
#S2 prepare input for MeQTL
cd $out_dir/$group
#bash /share/pub/xingsl/shilai/pipeline/tools/DaPars2/3aQTL-pipe/src/prepare_inputs_for_3aeQTL_mapping_forS1960.sh -g $out_dir/$group/vcf_list.txt -p $out_dir/$group/phenotype/$group\.paired.cpm.changeID.txt -c /APA/MeQTL/$group/phenotype/$group\.subpop -s $out_dir/$group/phenotype/$group\.lst -m 0.05 -n 5

/APA/software/miniforge3/envs/R3.6.3/bin/Rscript /share/pub/xingsl/shilai/pipeline/tools/DaPars2/3aQTL-pipe/src/curate_pheno_geno_covariates_forS1960_aeQTL.R -p $out_dir/$group/phenotype/$group\.paired.cpm.changeID.txt -c /APA/MeQTL/$group/phenotype/$group\.subpop -n 5 --cpm_threshold 5 --sample_pct 0.3

#S3 3aQTL_mapping
cd $out_dir/$group
/APA/software/miniforge3/envs/R3.6.3/bin/Rscript  /share/pub/xingsl/shilai/pipeline/tools/DaPars2/3aQTL-pipe/src/run_3aQTL_mapping.R -q 1e-2 -Q 1e-6
#S4 QTL plot visualize_sig
#/APA/software/miniforge3/envs/R3.6.3/bin/Rscript /share/pub/xingsl/shilai/pipeline/tools/DaPars2/3aQTL-pipe/src/QTL_plot.R -s "Chr10_132012791_C_T" -g "Zm00001d025866_T001" 
#S5 pre_fine-mapping, 用0.05/转录本数目
bash /share/pub/xingsl/shilai/pipeline/tools/DaPars2/3aQTL-pipe/src/prepare_inputs_for_finemapping.modi.sh -q 0.05
#S6 fine_mapping
#bash /share/pub/xingsl/shilai/pipeline/tools/DaPars2/3aQTL-pipe/src/run_fine_mapping.sh -t 12
bash /share/pub/xingsl/shilai/pipeline/tools/DaPars2/3aQTL-pipe/src/run_fine_mapping_modi.sh -t 12
#S7 merge_finemap_results
/APA/software/miniforge3/envs/R3.6.3/bin/Rscript /share/pub/xingsl/shilai/pipeline/tools/DaPars2/3aQTL-pipe/src/merge_finemap_results.R
