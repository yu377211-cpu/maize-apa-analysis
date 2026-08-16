#!usr/bin/python3
import argparse
from collections import defaultdict

# 创建解析器
parser = argparse.ArgumentParser(description="Process some integers.")
# 添加参数
parser.add_argument("input_file", type=str, help="Path to the input file")
parser.add_argument("output_file", type=str, help="Path to the output file")
# 解析参数
args = parser.parse_args()

#F = open("/share/pub/xingsl/shilai/project/APA/PAC_3seq_hq/region_dis/SRP115041_ab22.intersect.3UTR.annoTranscript.region.ratio.combine","r")
with open(args.input_file, "r") as F:
  lines = F.readlines()

sample_list = lines[0].strip().split("\t")[4:]

out_dict = defaultdict(dict)
for i in range(len(sample_list)):
  flag_dict = defaultdict(lambda:([],[],[]))
  for l in lines[1:]:
    transcript = l.split("\t")[0].split("|")[0]
    flag =  l.split("\t")[2]  
    #s = int(l.split("\t")[i+4])
    s = float(l.split("\t")[i+4])
    if flag == "S":
      flag_dict[transcript][0].append(s)
    elif flag == "L": 
      flag_dict[transcript][1].append(s)
    elif flag == "M":
      flag_dict[transcript][2].append(s)
  for f in flag_dict:  
    if sum(flag_dict[f][0]) + sum(flag_dict[f][1]) + sum(flag_dict[f][2]) == 0:
      r = "NA" 
    elif sum(flag_dict[f][1]) == 0:
      r = 0
    else :
      r = sum(flag_dict[f][1])/(sum(flag_dict[f][0]) + sum(flag_dict[f][1]) + sum(flag_dict[f][2]))     
    out_dict[f][sample_list[i]] = r


#o = open("/share/pub/xingsl/shilai/project/APA/PAC_3seq_hq/region_dis/SRP115041_ab22.intersect.3UTR.PDUI","w") 
with open(args.output_file, "w") as o:
  l = "\t".join(sample_list)
  o.write(f"Transcript\t{l}\n") 
  for i in out_dict:
    l = []
    for s in sample_list:
      l.append(str(out_dict[i][s]))
    L = "\t".join(l)
    o.write(f"{i}\t{L}\n")
#o.close()
