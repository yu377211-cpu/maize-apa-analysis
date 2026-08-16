cat SRP115041_ab22.intersect.3UTR SRP115041_ab22.intersect.5UTR SRP115041_ab22.intersect.CDS SRP115041_ab22.intersect.extend SRP115041_ab22.intersect.intron |awk -F"\t" '{if($7 != "."){print}}' >SRP115041_ab22.intersect.allGeneType.txt
cut -f5,10 SRP115041_ab22.intersect.allGeneType.txt|awk -F"_" '{print $1}'|sort -u >PAC_Gene.pair.txt
awk -F'\t' '{ if (arr[$2]) arr[$2] = arr[$2] "\t" $1; else arr[$2] = $1 } END { for (key in arr) print key "\t" arr[key] }' PAC_Gene.pair.txt > PAC_Gene.merged.txt
awk -F"\t" '{if($7 != "."){print}}' SRP115041_ab22.intersect.3UTR|cut -f5|sort -u >SRP115041_ab22.intersect.3UTR.PAC.lst
awk -F"\t" '{if($7 != "."){print}}' SRP115041_ab22.intersect.3UTR|cut -f5,10|awk -F"_" '{print $1}'|sort -u >SRP115041_ab22.intersect.3UTR.annoGene
perl /APA/PAC_3seq_diffana/bin/tiqu.pl SRP115041_ab22.intersect.3UTR.PAC.lst /APA/PAC_3seq_hq/region_dis/SRP115041_ab22.intersect.3UTR.tag SRP115041_ab22.intersect.3UTR.tag
python3 pac2gene.py -p /APA/PAC_3seq_hq/umap.Chr.PAC.PATcount.SampleExpressRatio_ab22.matrix -g PAC_Gene.merged.txt -o SRP115041_ab22.intersect.gene_level.tag
perl /APA/PAC_3seq_hq/tiqu_PAC.pl SRP115041_ab22.intersect.3UTR.PAC.lst PAC_Gene.pair.txt PAC_Gene.pair3UTR.txt
cut -f2 PAC_Gene.pair3UTR.txt|sort -u >PAC_Gene.3UTR.lst
grep -f PAC_Gene.3UTR.lst SRP115041_ab22.intersect.gene_level.tag >SRP115041_ab22.intersect.3UTR.gene_level.tag
