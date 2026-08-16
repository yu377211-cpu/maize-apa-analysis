dir=/APA/AFreview/PAC_3seq_diff_usage/same_subpop
bin=/APA/PAC_3seq_diff_usage/same_subpop/bin
subset=ln
#lst=" GRootvsGShoot GRootvsKern GRootvsL3Base GRootvsL3Tip GRootvsLMAD GRootvsLMAN GShootvsKern GShootvsL3Base GShootvsL3Tip GShootvsLMAD GShootvsLMAN KernvsL3Base KernvsL3Tip KernvsLMAD KernvsLMAN L3BasevsL3Tip L3BasevsLMAD L3BasevsLMAN L3TipvsLMAD L3TipvsLMAN LMADvsLMAN "
lst=" L3BasevsKern L3TipvsKern LMADvsKern LMANvsKern "

cd $dir
mkdir $subset
cd $dir/$subset


#sed -i s/sample/transcript/g $dir/$subset\.hq.3UTR.PDUI

export R_LIBS="/share/apps/R-3.6.3/lib64/R/library/":$R_LIBS
export R_LIBS="/APA/software/miniforge3/envs/R3.6.3/lib/R/library/":$R_LIBS
for i in ${lst}
do
	ln -s /APA/AFreview/PAC_3seq_diffana/PAC_readcount/same_subpop/diff/${subset}/$i\.condition .
	/APA/software/miniforge3/envs/R3.6.3/bin/Rscript $bin/limma.without_heatmap_block.R $dir/${subset}.hq.3UTR.PDUI ${i}.condition ${subset}.${i}
	/APA/software/miniforge3/envs/R3.6.3/bin/Rscript $bin/vocal.R  ${subset}.${i}.full.txt 0.05 0.2
	/APA/software/miniforge3/envs/R3.6.3/bin/Rscript $bin/heatmap.R $dir/${subset}.hq.3UTR.PDUI  ${subset}.${i}.sig.txt ${i}.condition ${subset}.${i}
done

#显著差异的转录本(adj.P.Val < 0.05 & |logFC| > 0.2)
