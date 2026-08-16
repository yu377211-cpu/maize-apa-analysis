#!/bin/bash

for i in GRoot GShoot GShoot Kern L3Base L3Tip LMAD LMAN
do
	cd /APA/AFreview/3aQTL/${i}/FineMapping
	Rscript /APA/AFreview/3aQTL/bin/summary_susie.R susieR_res.all_genes.txt susie_summary.txt
done

#paste -d "\t" /APA/AFreview/3aQTL/GRoot/FineMapping/susie_summary.txt /APA/AFreview/3aQTL/GShoot/FineMapping/susie_summary.txt /APA/AFreview/3aQTL/Kern/FineMapping/susie_summary.txt /APA/AFreview/3aQTL/L3Base/FineMapping/susie_summary.txt /APA/AFreview/3aQTL/L3Tip/FineMapping/susie_summary.txt /APA/AFreview/3aQTL/LMAD/FineMapping/susie_summary.txt /APA/AFreview/3aQTL/LMAN/FineMapping/susie_summary.txt >tmp.xls
