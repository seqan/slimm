#!/usr/bin/env python
import sys
import pandas as pd
import urllib
Name = []
Taxid = []
NoOfReads = []
RelativeAbundance = []
Contributers = []
Coverages = []

trim_left = 15
samples = {}
column_names = ["level", "taxid", "name", "linage"]
sample_count = 0

for file_path in sys.argv[1:]:
    f_name_start = file_path.rfind('/')
    dot_pos = file_path.rfind('.')
    samples[file_path[f_name_start+1:dot_pos]] = sample_count
    sample_count += 1
    column_names.append(file_path[f_name_start+1:dot_pos])

all_taxids = {}
# collect taxids from all files  
for file_path in sys.argv[1:]:
    inpf = open(file_path, 'r')
    first = True
    for line in inpf:
        if first :
            first = False
            continue
        values = line.split("\t")
        all_taxids[values[2]] = [values[0], values[1], values[3], values[2]]
    inpf.close()

# fill all the abundunce values with zeros
for  t_id in all_taxids:
    all_taxids[t_id] += sample_count*["0.0"]

# # fill all the actual abundunce
sample_no = 0
for file_path in sys.argv[1:]:
    inpf = open(file_path, 'r')
    first = True
    for line in inpf:
        if first :
            first = False
            continue
        values = line.split("\t")
        all_taxids[values[2]][4+sample_no] = values[4]
    sample_no += 1
    inpf.close()

sort_columns = ["level"] + column_names[3:]
sort_direct = len(sort_columns) * [False]
merged_out = pd.DataFrame.from_dict(all_taxids, orient='index')
merged_out.columns = column_names
merged_out.sort_values(sort_columns, ascending=sort_direct, inplace=True)

merged_out.to_csv(path_or_buf="merged_profile.tsv", sep='\t', index=False)
