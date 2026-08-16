for i in GRoot GShoot Kern L3Base L3Tip LMAD LMAN
do
        Rscript ../bin/run_PAC_coloc_susie.R --ae /APA/AFreview/MeQTL/${i}/FineMapping/lead_snps_list.txt --a /APA/AFreview/3aQTL/${i}/FineMapping/lead_snps_list.txt --map /APA/AFreview/PAC_3seq_hq/region_dis_addStrand/SRP115041_ab22.intersect.3UTR.annoTranscript --out ${i}_coloc
done

python3 ../bin/two_direction_bar.py -i lead_snp_pairs.stat -o lead_snp_pairs.proportion_plot
sh combine_detail.sh
