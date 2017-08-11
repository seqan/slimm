#!/usr/bin/env python
import argparse
import glob
import os
import re
import gzip

parser = argparse.ArgumentParser(description =
''' Download reference genomes of microorganisms
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('workdir', type=str,
                    help = 'The path of working directory where (intermidiate) results will be saved')

args = parser.parse_args()

working_dir = args.workdir
groups = ""
db_choice = ""


extension = ".fna.gz"
genomes_dir = ""
for file in os.listdir(working_dir):
    if "genomes" in file and not "genomes_to_download" in file:
        genomes_dir = file
    elif "genomes_to_download" in file:
        groups = file
    elif "assembly_summary_" in file:
        db_choice = file


groups = groups.replace("_genomes_to_download", "")
groups = groups.replace(".txt", "")

db_choice = db_choice.replace("assembly_summary_", "")
db_choice = db_choice.replace(".txt", "")

genomes_dir_path = working_dir + "/" + genomes_dir
combined_genomes_file = groups + "_" + db_choice + "_combined.fna"
combined_genomes_file_path = working_dir + "/" + combined_genomes_file


outf = open(combined_genomes_file_path, 'w')
files =  sorted(glob.glob(genomes_dir_path + "/*" + extension))
for fasta_file in files :
    taxon = fasta_file.replace(genomes_dir_path, "").replace(extension, "").replace("/", "")
    inpf = gzip.open(fasta_file, 'rb')
    current_id = ""
    current_seq = ""
    count = 0
    plasmid_count = 0
    line_len = 0
    is_plasmid = False
    for line in inpf:
        if re.search("\>", line):
            is_plasmid = "plasmid" in line.lower()
            if is_plasmid:
                plasmid_count += 1
                continue
            if count == 0:
                splited_line = line.split()
                current_id = splited_line[0] +"|kraken:taxid|" + taxon + " " +" ".join(splited_line[1:]) + "\n"
                outf.write(current_id)
            else:
                current_id = (line_len * "N") +"\n"
                outf.write(current_id)
            current_seq = ""
            count += 1
        elif not is_plasmid :
            outf.write(line)
            if line_len == 0:
                line_len = len(line) -1
    print taxon + extension ,".\t ", count , " seqs\t", plasmid_count ,"plasmids\tall seqs writeen delimated by aline of N's. all_plasmids are ignored"
    inpf.close()
outf.close()
print "merged file written to" , combined_genomes_file_path
