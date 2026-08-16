#!/usr/bin/python
import argparse
#按列的id提取列，且删除一行全为NA的情况
parser = argparse.ArgumentParser(description="find file row")
parser.add_argument("-k", "--keyfile", help="keyword file", required=True)
parser.add_argument("-f", "--filename", help="search file", required=True)
parser.add_argument("-o", "--outfile", help="output file", required=True)
argv = parser.parse_args()

key_file = argv.keyfile.strip()
search_file = argv.filename.strip()
o_file = argv.outfile.strip()

fn = open(key_file, "r", encoding="utf-8")

key = []
for line in fn:
    key.append(line.strip())
fn.close()
# 保留第一列
key.insert(0, "Transcript")

fp = open(search_file, "r")
op = open(o_file, "w")

title = fp.readline().strip().split("\t")
op.writelines("\t".join(key) + "\n")

for line in fp:
    line_list = line.strip().split("\t")
    line_dict = dict(zip(title, line_list))
    write_list = []
    for item in key:
        if line_dict.get(item):
            write_list.append(line_dict.get(item))
    #if list(filter(lambda x: float(x) > 0, write_list[1:])):
    if not all(item == "NA" for item in write_list):
        op.writelines("\t".join(write_list) + "\n")

fp.close()
op.close()
