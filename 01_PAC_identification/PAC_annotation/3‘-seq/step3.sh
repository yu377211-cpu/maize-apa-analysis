Rscript cal_TMM_CPM.R umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.matrix
perl /APA/PAC_3seq_hq/tiqu_PAC.pl SRP115041_ab22.intersect.3UTR.PAC.lst umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.matrix_TMM_CPM.txt umap.Chr.ab22.matrix_TMM_CPM.3UTR.txt
perl /APA/PAC_3seq_hq/tiqu_PAC.pl SRP115041_ab22.intersect.3UTR.PAC.lst umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.matrix_TMM_log2CPM.txt umap.Chr.ab22.matrix_TMM_log2CPM.3UTR.txt
Rscript plot_PCA_forCPM.r -e umap.Chr.ab22.matrix_TMM_CPM.3UTR.txt -m /APA/PAC_3seq_diffana/S1960.group.txt -o CPM_3UTR_PCA
Rscript plot_CPM_distribution.R umap.Chr.ab22.matrix_TMM_CPM.3UTR.txt umap.Chr.ab22.matrix_TMM_CPM.3UTR
