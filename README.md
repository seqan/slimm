SLIMM - Species Level Identification of Microbes from Metagenomes - 
=================================================================

a taxonomic profiling tool that investigates which microorganisms are present in a sequenced sample. SLIMM requires a BAM/SAM alignment file as an input. One can use a read mapper of choice to map raw reads obtained from a sequencing machine to obtain the BAM/SAM file required as input for SLIMM. 

    slimm [OPTIONS] -m $SLIMM_DB_PATH $SAM_FILE_PATH
    Try 'slimm --help' for more information.

VERSION

    * SLIMM version: 0.2.1
    * Last update: January 2017
    

### Downloads:

#### Executables
Pre-built executables for Linux and Mac are made available at the [releases page]( https://github.com/seqan/slimm/releases/latest).

#### Source code
You can build SLIMM from its source. Instruction on how to build from source can be found at the [slimm wiki] (https://github.com/seqan/slimm/wiki) 

#### Cite us

If you use SLIMM in your work-flows, don't forget to cite us.

Dadi TH, Renard BY, Wieler LH, Semmler T, Reinert K. (2017) SLIMM: species level identification of microorganisms from metagenomes. PeerJ 5:e3138 https://doi.org/10.7717/peerj.3138
