#!/bin/bash

lst=( ln mixed nss popcorn ss sweet ts )
bin=/APA/PAC_3seq_diff_usage/same_subpop/bin

for i in ${lst[@]}
do
	python3 $bin/find_row_pdui.py -k /APA/PAC_3seq_diffana/PAC_tpm_ab22/same_subpop/${i}.lst -f /APA/AFreview/PAC_3seq_hq/region_dis_addStrand/PDUI/SRP115041_ab22.intersect.3UTR.annoTranscript_ratio.tag.PDUI_rmdup -o $i\.hq.3UTR.PDUI
	sed -i s/Transcript/sample/g $i\.hq.3UTR.PDUI
	awk 'BEGIN{sum = 0} {for(i = 2; i <= NF; i++) {if($i == "NA") sum += 1} {print $1"\t"sum; sum = 0}}' $i\.hq.3UTR.PDUI|cut -f2|sort|uniq -c|sort -k2n|awk '{print $2"\t"$1}' >$i\.hq.3UTR.PDUI.SampleUsage.sta
	#/share/apps/R-4.2.1/bin/Rscript /APA/PAC_3seq_diffana/bin/num_sta.ratio.R $i\.hq.3UTR.PDUI.SampleUsage.sta $i\.3UTR.SampleUsage.ratio.pdf $i
done
