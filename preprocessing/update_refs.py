#!/usr/bin/env python
from helper_methods import *

parser = argparse.ArgumentParser(description = 
''' Download reference genomes of microorganisms  
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-wd', '--workdir', type=str, 
                    help = 'The path of working directory where (intermidiate) results will be saved')
parser.add_argument('-g', '--groups',  type=str, default = "AB", 
                    help = '''Which group of microbes to consider any combination of the letters [A], [B] and [V] 
                    where B =  Bacteria, A = Archaea and V = Viruses and Viroids (default: AB)''')
parser.add_argument('-s', '--sp', dest='species_lv', action='store_true',
                    help = 'download one reference per species.')
parser.add_argument('-t', '--threads',  type=int, choices=xrange(1, 11), default = 4, 
                    help = 'number of threads for downloading in parallel in the range 1..10 (default: 4)')

args = parser.parse_args()

parallel = args.threads
working_dir = args.workdir
groups = args.groups
only_species = args.species_lv
db_choice = "refseq"
old_date_string = ""

exists_genomes_to_download = False
exists_genomes = False
exists_assembly_summary = False
exists_slimmDB = False
exists_taxcat = False
exists_taxdump = False

old_genomes_to_download_path = ""
old_slimmDB_path = ""
old_genomes_dir = ""

for file in os.listdir(working_dir):
    if file.endswith(".old") :
        continue
    if "genomes_to_download" in file:
        exists_genomes_to_download = True
        old_genomes_to_download_path = working_dir + "/.old/" + file
    elif "genomes" in file:
        exists_genomes = True
        old_genomes_dir =  working_dir + "/.old/" + file
    elif "assembly_summary" in file:
        exists_assembly_summary = True
        if "assembly_summary_genbank_" in file :
            db_choice = "genbank"
            old_date_string = file[25:33]
        else:
            old_date_string = file[24:32]
    elif "slimmDB_" in file:
        exists_slimmDB = True
        old_slimmDB_path = working_dir + "/.old/" + file
    elif "taxcat_" in file:
        exists_taxcat = True
    elif "taxdump_" in file:
        exists_taxdump = True

download_update = False
if exists_genomes_to_download and exists_genomes and exists_assembly_summary and exists_slimmDB and exists_taxcat and exists_taxdump :
    existing_date = datetime.datetime.strptime(old_date_string, '%d%m%Y')
    diff = datetime.datetime.now() - existing_date
    download_update = query_yes_no("The database is [" + str(diff.days) + "] days old. Do you realy want to update it?")
    if not download_update:
        sys.exit(0)
    if not os.path.isdir(working_dir + "/.old/"):
        os.makedirs(working_dir + "/.old/")
    for file in os.listdir(working_dir):   
        if file.endswith(".old") :
            continue
        os.rename(working_dir + "/" + file, working_dir + "/.old/" + file)
else:
    print "Error! corrupted database!" 
    sys.exit(1)

try:
    today_string = "25082016"
    # today_string = (datetime.datetime.now()).strftime("%d%m%Y")
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
    taxaid2sp_path=nodes_path.replace("nodes.dmp", "taxaid2sp_" + groups + ".dmp")


    assembly_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
    if db_choice == "genbank":
        assembly_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"

    assembly_summary_file = working_dir + "/assembly_summary_" + db_choice + "_" + today_string + ".txt"

    genomes_to_download_path = working_dir + "/" + groups + "_genomes_to_download.txt"
    ncbi_get_script_path =  working_dir + "/" + groups + "_genomes_ncbi_download.sh"

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
        if l[0] in groups :
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
                fname = path.replace("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/", "")
                complete_path = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/" + fname + "/" + fname + "_genomic.fna.gz"
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


    old_urls = {}
    new_urls = {}
    updated_urls = {}
    added_urls = {}

    inpf = open(old_genomes_to_download_path, 'r')

    # DELETE FOR TESTING ONLY 
    # count = 0
    # for line in inpf:
    #     if count == 96:
    #         break
    #     count += 1
    for line in inpf:
        l = line.split('\t')
        old_urls[l[0]] = l[4].strip()
    inpf.close()

    inpf = open(genomes_to_download_path, 'r')

    # DELETE FOR TESTING ONLY 
    # count = 0
    # for line in inpf:
    #     if count == 100:
    #         break
    #     count += 1
    for line in inpf:
        l = line.split('\t')
        new_urls[l[0]] = l[4].strip()
    inpf.close()

    for taxid in new_urls:
        if taxid in old_urls:
            if new_urls[taxid] == old_urls[taxid] :
                old_file_path = old_genomes_dir + "/" + str(taxid) + ".fna.gz"
                new_file_path = genomes_dir + "/" + str(taxid) + ".fna.gz"
                shutil.copyfile(old_file_path, new_file_path)
            else :
                updated_urls[taxid] =  new_urls[taxid]
        else:
            added_urls[taxid] =  new_urls[taxid]
    print (str(len(updated_urls)) + " references will be updated and " + str(len(added_urls)) + " will be added.")

    ########################################################################################################################
    # Reduce the nodes, and names lookup files so that they contain nodes and names of groups of interest (BAV in this case)
    # and make sure that all the nodes from the old database are included
    # B =  Bacteria
    # A = Archaea
    # V = Viruses and Viroids
    ########################################################################################################################

    taxid_parent =  {} 
    taxid_rank = {}
    names = {}
    # load old nodes to taxid_parent and taxid_rank
    inpf = open(old_slimmDB_path + "/nodes.dmp" , 'r')
    for line in inpf:
        l = line.split('\t')
        taxid = int(l[0])
        parent = int(l[1])
        rank = l[2]
        taxid_parent[taxid] = parent
        taxid_rank[taxid] = rank
    inpf.close()
    # load old names to dict names
    inpf = open(old_slimmDB_path + "/names.dmp" , 'r')
    for line in inpf:
        l = line.split("\t")
        names[int(l[0])] = l[1]
    inpf.close()

    print "Reducing and updating nodes to interest groups i.e. [" + groups + "]"

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
    print "results (nodes mapping ) are written to " + reduced_nodes_path
    print "results (name mapping ) are written to " + reduced_names_path

    print "Updating reference genomes. This might take a while! ..."
    download_queue = Queue.Queue()
    for ti in updated_urls:
        download_item = [ti, updated_urls[ti]]
        download_queue.put(download_item) # produce

    for ti in added_urls:
        download_item = [ti, added_urls[ti]]
        download_queue.put(download_item) # produce

    total_to_download = len(updated_urls) + len(added_urls)
    while not download_queue.empty():
        if threading.activeCount() <= parallel:
            threading.Thread(target=download_one, args=(download_queue, genomes_dir)).start()
            sys.stdout.write('\r')
            sys.stdout.write("%d of %d remaining ..." % (download_queue.qsize(), total_to_download))
            sys.stdout.flush()
except Exception, e:
    print e
    print "Update not complete! rolling back changes!"
    for file in os.listdir(working_dir):
        if file.endswith(".old") :
            continue
        elif os.path.isdir(working_dir + "/" + file):
            shutil.rmtree(working_dir + "/" + file)
        else :
            os.remove(working_dir + "/" + file)
    for file in os.listdir(working_dir + "/.old/"):
        shutil.move(working_dir + "/.old/" + file, working_dir + "/" + file)
    sys.exit(0)

# delete the old files
shutil.rmtree(working_dir + "/.old/")