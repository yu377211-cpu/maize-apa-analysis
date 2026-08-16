perl distance_ratio.pl /APA/AFreview/PAC_3seq_hq/region_dis_addStrand/SRP115041_ab22.intersect.3UTR SRP115041_ab22.intersect.3UTR_region.ratio.old 
perl distance_ratio_mulUTR.pl /APA/PAC_3seq_hq/region_dis/three_prime_UTR.bed /APA/AFreview/PAC_3seq_hq/region_dis_addStrand/SRP115041_ab22.intersect.3UTR SRP115041_ab22.intersect.3UTR_region.ratio
perl combine.pl ../SRP115041_ab22.intersect.3UTR.tag SRP115041_ab22.intersect.3UTR_region.ratio SRP115041_ab22.intersect.3UTR.annoTranscript_ratio.tag
python3 /APA/PAC_3seq_hq/region_dis/addM/calculate_PDUI_addM.py SRP115041_ab22.intersect.3UTR.annoTranscript_ratio.tag SRP115041_ab22.intersect.3UTR.annoTranscript_ratio.tag.PDUI
perl /APA/PAC_3seq_hq/region_dis/rmdup.pl SRP115041_ab22.intersect.3UTR.annoTranscript_ratio.tag.PDUI >SRP115041_ab22.intersect.3UTR.annoTranscript_ratio.tag.PDUI_rmdup
gzip -c SRP115041_ab22.intersect.3UTR.annoTranscript_ratio.tag.PDUI_rmdup >hq_3UTR_PDUI_3seq.tsv.gz 
Rscript PDUI_distribute_plot.r
