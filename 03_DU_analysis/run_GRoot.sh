dir=/APA/AFreview/PAC_3seq_diff_usage/same_tissue
bin=/APA/PAC_3seq_diff_usage/same_subpop/bin
subset=GRoot
lst=" mixedvsnss mixedvspopcorn mixedvsss mixedvssweet mixedvsts nssvspopcorn nssvsss nssvssweet nssvsts popcornvsss popcornvssweet popcornvsts ssvssweet ssvsts sweetvsts lnvsmixed lnvsnss lnvspopcorn lnvsss lnvssweet lnvsts "
cd $dir
mkdir -p $subset
cd $dir/$subset


#awk 'END{print NR-1}' $dir/../$subset\.hq.3UTR.PDUI
#20327
#sed -i s/sample/transcript/g $dir/../$subset\.hq.3UTR.TPM

export R_LIBS="/share/apps/R-3.6.3/lib64/R/library/":$R_LIBS
export R_LIBS="/APA/software/miniforge3/envs/R3.6.3/lib/R/library/":$R_LIBS
for i in ${lst}
do
	ln -s /APA/AFreview/PAC_3seq_diffana/PAC_readcount/same_tissue/diff//${subset}/${i}.condition .
	/APA/software/miniforge3/envs/R3.6.3/bin/Rscript $bin/limma.without_heatmap.R $dir/${subset}.hq.3UTR.PDUI ${i}.condition ${subset}.${i}
	/APA/software/miniforge3/envs/R3.6.3/bin/Rscript $bin/vocal.R  ${subset}.${i}.full.txt 0.05 0.2
	/APA/software/miniforge3/envs/R3.6.3/bin/Rscript $bin/heatmap.R $dir/${subset}.hq.3UTR.PDUI  ${subset}.${i}.sig.txt ${i}.condition ${subset}.${i}
done

#显著差异的转录本(adj.P.Val < 0.05 & |logFC| > 0.2)
