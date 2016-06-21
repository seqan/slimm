SLIMM - Species Level Identification of Microbes from Metagenomes - Investigates which microbes are present from a BAM/SAM alignment file.
=================================================================
    slimm [OPTIONS] "IN" "OUT"
    Try 'slimm --help' for more information.

VERSION

    * SLIMM version: 0.1
    * Last update: June 2016
    
### How do I get set up? ###

In order to run SLIMM you need the following files which are made available at http://ftp.mi.fu-berlin.de/pub/dadi/slimm/.

<ol type="1">
	<li> a bowtie2 index of bacterial and archeal genomes. Here you have two options:
		<ol type="a">	
			<li> a reference genome database of 4538 complete bacterial and archial genomes that covers around 2500 different species (AB_complete.tar.gz) </li>
			<li> a reference genome database of 13192 complete and draft bacterial and archial genomes. Here 1 species is represented by 1 genome. (AB_species.tar.gz)</li>
		</ol>
	</li>
	<li> a SLIMM database based on your choice of the reference genome database:
		<ol type="a">
	  		<li> slimmDB-4538.tar.gz </li>
	  		<li> slimmDB-13192.tar.gz </li>
		</ol>
	</li>
	<li> the binary executable (slimm) </li>
	<li> bowtie2 mapper </li>
</ol>

### Example ###

We assume you have the following folder/file structure in your working directory:
```
Working Directory
  │
  ├── slimm (slimm binary executable)
  │
  ├── bowtie2_indices (indexed reference genomes)
  ├    ├── AB_complete
  ├    ├    ├── AB_complete.1.bt2l
  ├    ├    ├── AB_complete.2.bt2l
  ├    ├    ├── ...
  ├    ├── AB_species
  ├    ├    ├── AB_species.1.bt2l
  ├    ├    ├── AB_species.2.bt2l
  ├    ├    ├── ...
  ├── slimmDB-13192 (SLIMM taxonomic database to be used with AB_species)
  │
  ├── slimmDB-4538 (SLIMM taxonomic database to be used with AB_complete)
  │
  ├── alignment_files (alignment files will be stored here)
  │
  ├── slimm_reports (slimm taxonomic reports will be stored here)
  │
  ├── mg_reads (metagenomic sequencing reads) 
  ├    ├── SRR1748536_1.fastq
  ├    ├── SRR1748536_2.fastq
```  

Use bowtie2 to map the metagenomic reads against reference genomes and produce alignment files.

	bowtie2 -x ./AB_species/AB_species -1 ./mg_reads/SRR1748536_1.fastq -2 ./mg_reads/SRR1748536_1.fastq -q --no-unal --mm -p 10 -k 60 -S ./alignment_files/SRR1748536.sam

[recomended] Alternatively we can directly pipe the output to samtools to save it as BAM format. This is faster than just writing into sam files since it avoids a lot of io operations. And SLIMM performs sligtly better on BAM format.
	
	bowtie2 -x ./AB_species/AB_species -1 ./mg_reads/SRR1748536_1.fastq -2 ./mg_reads/SRR1748536_1.fastq -q --no-unal --mm -p 10 -k 60 2>./alignment_files/SRR1748536.txt | samtools view -bS - > ./alignment_files/SRR1748536.bam;

Run SLIMM on the output of the read mapper (SAM/BAM files)
	
	slimm -m ./slimmDB-13192 ./alignment_files/SRR1748536.bam -o slimm_reports/

you will find two reports under the directory slimm_reports. 

	* SRR1748536.tsv contains raw information about all individual genomes in the database and 
	* SRR1748536_sp_reported.tsv contains the species that SLIMM reported to be present in the sample.

You can also tell SLIMM to report at lower level of the the taxonomic tree genus family ... (see slimm --help for more details)

<!---

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines --->

### Who do I talk to? ###

* Temegsen H. Dadi 

    email: temesgen.dadi[at]fu-berlin.de

* The SeqAn team 

    website: www.seqan.de
