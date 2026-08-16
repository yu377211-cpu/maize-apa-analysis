#lst=" GRootvsGShoot GRootvsKern GRootvsL3Base GRootvsL3Mid GRootvsL3Tip GRootvsLMAD GRootvsLMAN GShootvsKern GShootvsL3Base GShootvsL3Mid GShootvsL3Tip GShootvsLMAD GShootvsLMAN KernvsL3Base KernvsL3Mid KernvsL3Tip KernvsLMAD KernvsLMAN L3BasevsL3Mid L3BasevsL3Tip L3BasevsLMAD L3BasevsLMAN L3MidvsL3Tip L3MidvsLMAD L3MidvsLMAN L3TipvsLMAD L3TipvsLMAN LMADvsLMAN "
lst=" GRootvsGShoot GRootvsKern GRootvsL3Base GRootvsL3Tip GRootvsLMAD GRootvsLMAN GShootvsKern GShootvsL3Base GShootvsL3Tip GShootvsLMAD GShootvsLMAN L3BasevsKern L3TipvsKern LMADvsKern LMANvsKern L3BasevsL3Tip L3BasevsLMAD L3BasevsLMAN L3TipvsLMAD L3TipvsLMAN LMADvsLMAN "
group=" ln mixed nss popcorn ss sweet ts "

cd /APA/AFreview/PAC_3seq_diff_usage/same_subpop

for j in ${group}
do
	for i in ${lst}
	do
		awk -F"\t" '{if($2 >0){print $1}}' ${j}/${j}.${i}.sig.txt |grep -v "Transcript" >${j}/${j}.${i}.sig.lengthening
		awk -F"\t" '{if($2 <0){print $1}}' ${j}/${j}.${i}.sig.txt |grep -v "Transcript" >${j}/${j}.${i}.sig.shortening
	done
done
