lst=" mixedvsnss mixedvspopcorn mixedvsss mixedvssweet mixedvsts nssvspopcorn nssvsss nssvssweet nssvsts popcornvsss popcornvssweet popcornvsts ssvssweet ssvsts sweetvsts lnvsmixed lnvsnss lnvspopcorn lnvsss lnvssweet lnvsts "
group=" GRoot GShoot Kern L3Base L3Tip LMAD LMAN "

cd /APA/PAC_3seq_diff_usage/same_tissue

for j in ${group}
do
	for i in ${lst}
	do
		awk -F"\t" '{if($2 >0){print $1}}' ${j}/${j}.${i}.sig.txt |grep -v "Transcript" >${j}/${j}.${i}.sig.lengthening
		awk -F"\t" '{if($2 <0){print $1}}' ${j}/${j}.${i}.sig.txt |grep -v "Transcript" >${j}/${j}.${i}.sig.shortening
	done
done
