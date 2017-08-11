#!/usr/bin/env python
from helper_methods import *

parser = argparse.ArgumentParser(description =
''' Download reference genomes of microorganisms
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-wd', '--workdir', type=str, required=True,
                    help = 'The path of working directory where (intermediate) results will be saved')
parser.add_argument('-g', '--groups',  type=str, default = "AB",
                    help = '''Which group of microbes to consider any combination of the letters [A], [B] and [V]
    where B =  Bacteria, A = Archaea and V = Viruses and Viroids (default: AB)''')
parser.add_argument('-s', '--sp', dest='species_lv', action='store_true',
                    help = 'download one reference per species.')
parser.add_argument('-t', '--taxa_ids',  type=str, default = "",
                    help = '''comma separated list of taxonomic ids to be included (in addition to --groups) into
    the reference database database. This way you might even add the genome of Eukaryotes.
    e.g. the host genome''')
parser.add_argument('-tr', '--threads',  type=int, choices=xrange(1, 11), default = 4,
                    help = 'number of threads for downloading in parallel in the range 1..10 (default: 4)')
parser.add_argument('-d', '--database', type=str, choices = ['refseq', 'genbank'], default = 'refseq',
                    help = 'From which database references should be downloaded  (default: refseq)')
parser.add_argument('-ts', '--testing',  dest='testing', action='store_true',
                    help = 'This is a test run. Download only few genomes.')


args = parser.parse_args()

parallel = args.threads
working_dir = args.workdir
groups = args.groups
only_species = args.species_lv
db_choice = args.database
testing = args.testing
taxid_list = args.taxa_ids


##############################################################################
# For KNIME workflow only
##############################################################################
# parallel        = flow_variables['threads']
# working_dir     = flow_variables['workdir']
# groups          = flow_variables['groups']
# only_species    = flow_variables['species_lv']
# db_choice       = flow_variables['database']
# testing         = flow_variables['testing']

subset_taxids   = []
if len(taxid_list) >= 1:
    subset_taxids   = map(int, taxid_list.split(','))

groups_name     = groups
if len(groups_name) < 1:
    groups_name = "CUSTOM"
elif len(subset_taxids) >= 1:
    groups_name += "_CUSTOM"

if not os.path.isdir(working_dir):
    os.makedirs(working_dir)
else:
    empty = os.listdir(working_dir) == [];
    if not empty:
        print ("[ERROR!] Working directory [" + working_dir + "] should be empty!")
        sys.exit(0)

today_string = (datetime.datetime.now()).strftime("%d%m%Y")
genomes_dir = working_dir + "/genomes_" + today_string
if not os.path.isdir(genomes_dir):
    os.makedirs(genomes_dir)

slimmDB_dir = working_dir + "/slimmDB_" + today_string
if not os.path.isdir(slimmDB_dir):
    os.makedirs(slimmDB_dir)

##############################################################################
# Get the path of Names, nodes, catagories file from the extracted folder
##############################################################################

taxdmp_extract_dir = taxonomy_download("taxdump", working_dir, today_string)
taxcat_extract_dir = taxonomy_download("taxcat", working_dir, today_string)


names_path = taxdmp_extract_dir + "/names.dmp"
nodes_path = taxdmp_extract_dir + "/nodes.dmp"
catagories_path = taxcat_extract_dir + "/categories.dmp"
reduced_names_path = slimmDB_dir + "/names.dmp"
reduced_nodes_path = slimmDB_dir + "/nodes.dmp"
taxaid2sp_path=nodes_path.replace("nodes.dmp", "taxaid2sp_" + groups_name + ".dmp")


assembly_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
if db_choice == "genbank":
    assembly_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"

assembly_summary_file = working_dir + "/assembly_summary_" + db_choice + "_" + today_string + ".txt"

genomes_to_download_path = working_dir + "/" + groups_name + "_genomes_to_download.txt"
ncbi_get_script_path =  working_dir + "/" + groups_name + "_genomes_ncbi_download.sh"

print "Downloading assembly_summary file ..."
urllib.urlretrieve(assembly_summary_url, assembly_summary_file)


################################################################################################
# From the assembly summary file, identify, download and anotate the genomes to download.
# Unique genomes per taxon will be downloaded. Genomes are selected in the folowing order:
# 1. Complete genomes (longest first)
# 2. Chromosomes (longest first)
################################################################################################
taxid_col = 5
if only_species :
    taxid_col = 6

intial_taxids = {}
inpf = open(catagories_path, 'r')
for line in inpf:
    l = line.split('\t')
    if l[0] in groups or int(l[1]) in subset_taxids or  int(l[2]) in subset_taxids:
        intial_taxids[int(l[1])] = 1
        intial_taxids[int(l[2])] = 1
inpf.close()
taxid_genomes =  {}
inpf = open(assembly_summary_file, 'r')
firstLine = True
for line in inpf:
    if line[0] == "#":
        continue
    if ("representative genome" in line) or  ("reference genome" in line) or ("Complete Genome" in line) or ("Chromosome" in line) or ("Scaffold" in line) or ("Contig" in line):
    # if not firstLine:
        l = line.split('\t')
        taxid = int(l[taxid_col])
        path = l[19].replace("\n", "")
        if (taxid in intial_taxids and (path != 'na') and l[10] == "latest") :
            complete_path = path + path[path.rfind('/'):] + "_genomic.fna.gz"
            # fname = path.replace("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/", "")
            # complete_path = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/" + fname + "/" + fname + "_genomic.fna.gz"
            if taxid in taxid_genomes :
                taxid_genomes[taxid].append([l[11], l[13], complete_path, l[6], l[4]])
            else:
                taxid_genomes[taxid] = []
                taxid_genomes[taxid].append([l[11], l[13], complete_path, l[6], l[4]])
    firstLine = False
inpf.close()
outf = open(genomes_to_download_path, 'w')
outf_download = open(ncbi_get_script_path, 'w')
outf_download.write("#!/bin/bash\n\n")
downloader="~/local/bin/ncbi-get.sh "
downloader="ncbi-get.sh "
countxxx = 0
for taxid in taxid_genomes:
    found = False
    selected_op = []
    for options in taxid_genomes[taxid] :
        if (options[4] == "reference genome") :
            selected_op = options
            found = True
            break;
    if not found :
        for options in taxid_genomes[taxid] :
            if (options[4] == "representative genome") :
                selected_op = options
                found = True
                break;
    if not found :
        for options in taxid_genomes[taxid] :
            if (options[0] == "Complete Genome") :
                selected_op = options
                found = True
                break;
    if not found :
        for options in taxid_genomes[taxid] :
            if (options[0] == "Chromosome") :
                selected_op = options
                found = True
                break;
    if not found and (str(taxid) == options[3] or len(taxid_genomes[taxid]) == 1):
        for options in taxid_genomes[taxid] :
            if (options[0] == "Scaffold") :
                selected_op = options
                found = True
                break;
    if not found and (str(taxid) == options[3] or len(taxid_genomes[taxid]) == 1):
        for options in taxid_genomes[taxid] :
            if (options[0] == "Contig") :
                selected_op = options
                found = True
                break;
    if found:
        outf.write(str(taxid) + "\t" + selected_op[3] + "\t" + selected_op[0] +
                    "\t"+ selected_op[1] + "\t" + selected_op[2]  + "\n")
        outf_download.write(downloader + selected_op[2] + " " + genomes_dir + "/" + str(taxid) + ".fna.gz;")
outf.close()
outf_download.close()


########################################################################################################################
# Reduce the nodes, and names lookup files so that they contain nodes and names of groups of interest (BAV in this case)
# B =  Bacteria
# A = Archaea
# V = Viruses and Viroids
########################################################################################################################

print "Reducing nodes to interest groups i.e. [" + groups_name + "]"
taxid_parent =  {}
taxid_rank = {}
names = {}

inpf = open(names_path, 'r')
for line in inpf:
    if "scientific name" in line:
        l = line.split("\t|\t")
        names[int(l[0])] = l[1]
inpf.close()


inpf = open(nodes_path, 'r')
for line in inpf:
    l = line.split("\t|\t")
    taxid = int(l[0])
    parent = int(l[1])
    rank = l[2]
    taxid_parent[taxid] = parent
    taxid_rank[taxid] = rank
inpf.close()
print len(taxid_parent), "nodes found"
taxids_to_consider = {}
for taxid in intial_taxids:
    current_parent = taxid
    while(current_parent != 1):
        if current_parent in taxid_parent :
            taxids_to_consider[current_parent] = 1
            current_parent = taxid_parent[current_parent]
        else:
            print current_parent, " is not in the node file. may be, it is in deleted nodes!"
            break

outf = open(reduced_nodes_path, 'w')
outf2 = open(reduced_names_path, 'w')
for taxid in taxids_to_consider:
    outf.write(str(taxid) + "\t" + str(taxid_parent[taxid])+ "\t" + str(taxid_rank[taxid]) + "\n")
    outf2.write(str(taxid) + "\t" + names[taxid]  + "\n")
outf.close()
outf2.close()
print "results (nodes mapping ) are writen to " + reduced_nodes_path
print "results (name mapping ) are writen to " + reduced_names_path

print "Downloading reference genomes. This might take a while! ..."

download_queue = Queue.Queue()
inpf = open(genomes_to_download_path, 'r')

# DELETE FOR TESTING ONLY
# for line in inpf:
count = 0
for line in inpf:
    if count == 50 and testing:
        break
    count += 1
    l = (line.strip()).split("\t")
    download_item = [l[0], l[4]]
    download_queue.put(download_item) # produce
inpf.close()

total_to_download = download_queue.qsize()
while not download_queue.empty():
    if threading.activeCount() <= parallel:
        threading.Thread(target=download_one, args=(download_queue, genomes_dir)).start()
        sys.stdout.write('\r')
        sys.stdout.write("%d of %d remaining ..." % (download_queue.qsize(), total_to_download))
        sys.stdout.flush()