#!/usr/bin/python
import argparse

parser = argparse.ArgumentParser(description="Extract columns and filter rows with all NA values")
parser.add_argument("-k", "--keyfile", help="File containing column names to extract (with 'sample' header)", required=True)
parser.add_argument("-f", "--filename", help="Input file to process", required=True)
parser.add_argument("-o", "--outfile", help="Output file", required=True)
argv = parser.parse_args()

key_file = argv.keyfile.strip()
search_file = argv.filename.strip()
o_file = argv.outfile.strip()

# Read column names to extract (skip header line with 'sample')
with open(key_file, "r", encoding="utf-8") as fn:
    # Read all lines and skip the first line if it's 'sample'
    lines = fn.readlines()
    if lines and lines[0].strip().lower() == "sample":
        key = [line.strip() for line in lines[1:]]  # Skip header
    else:
        key = [line.strip() for line in lines]  # No header to skip
    
# Always keep the first column (Transcript)
if "Transcript" not in key:
    key.insert(0, "Transcript")

# Process input file
with open(search_file, "r") as fp, open(o_file, "w") as op:
    # Read and write header
    title = fp.readline().strip().split("\t")
    op.write("\t".join(key) + "\n")
    
    # Process each line
    for line in fp:
        line_list = line.strip().split("\t")
        line_dict = dict(zip(title, line_list))
        
        # Extract requested columns
        write_list = []
        for item in key:
            write_list.append(line_dict.get(item, "NA"))  # Use "NA" if column not found
        
        # Check if there's at least one non-NA value (excluding the first column)
        has_valid_data = any(val != "NA" for val in write_list[1:])
        
        if has_valid_data:
            op.write("\t".join(write_list) + "\n")
