dir=/APA/AFreview/PAC_3seq_diffana/PAC_readcount/same_subpop/diff
bin=/APA/PAC_3seq_diffana/bin
subset=ln
#lst=" GRootvsGShoot GRootvsKern GRootvsL3Base GRootvsL3Tip GRootvsLMAD GRootvsLMAN GShootvsKern GShootvsL3Base GShootvsL3Tip GShootvsLMAD GShootvsLMAN L3BasevsKern L3TipvsKern LMADvsKern LMANvsKern L3BasevsL3Tip L3BasevsLMAD L3BasevsLMAN L3TipvsLMAD L3TipvsLMAN LMADvsLMAN "
#lst=" L3BasevsKern L3TipvsKern LMADvsKern LMANvsKern "
lst=" LMANvsKern "

cd $dir
grep -f ../${subset}.lst /APA/PAC_3seq_diffana/S1960.tissue_genotype.txt >${subset}.group
mkdir -p $subset
cd $dir/$subset

for i in ${lst}
do
	python3  $bin/get_group_condition.covariate.py -i ${dir}/$subset\.group -c $i -o $i\.condition
done


export R_LIBS="/share/apps/R-3.6.3/lib64/R/library/":$R_LIBS
export R_LIBS="/APA/software/miniforge3/envs/R3.6.3/lib/R/library/":$R_LIBS
for i in ${lst}
do
	/APA/software/miniforge3/envs/R3.6.3/bin/Rscript $bin/DESeq2_covariate.R $dir/../$subset\.hq.3UTR.tag $i\.condition $subset\.$i\.heatmap.png "~ genotype + tissue"
done
