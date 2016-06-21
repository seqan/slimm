SLIMM - Species Level Identification of Microbes from Metagenomes - Investigates which microbes are present from a BAM/SAM alignment file.
=================================================================
    slimm [OPTIONS] "IN" "OUT"
    Try 'slimm --help' for more information.

VERSION

    * SLIMM version: 0.1
    * Last update: June 2016
    
### How do I get set up? ###

In order to run SLIMM you need the following files which are made available at http://ftp.mi.fu-berlin.de/pub/dadi/slimm/.

1. a bowtie2 index of bacterial and archeal genomes. Here you have two options:

<ol type="a">
  <li> a reference genome database of 4538 complete bacterial and archial genomes that covers around 2500 different species (AB_complete.tar.gz) </li>
  <li> a reference genome database of 13192 complete and draft bacterial and archial genomes. Here 1 species is represented by 1 genome. (AB_species.tar.gz)</li>
</ol>

2. a SLIMM database based on your choice of the reference genome database
<ol type="a">
  <li> slimmDB-4538.tar.gz </li>
  <li> slimmDB-13192.tar.gz </li>
</ol>
	
3. the binary executable (slimm)
4. bowtie2 mapper

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
