// ==========================================================================
//    SLIMM - Species Level Identification of Microbes from Metagenomes.
// ==========================================================================
// Copyright (c) 2014-2017, Temesgen H. Dadi, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Temesgen H. Dadi or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL TEMESGEN H. DADI OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
// ==========================================================================
#include "slimm.h"
using namespace seqan;

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("slimm");
    setDateAndVersion(parser);
    setDescription(parser);
    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIIN\\fP\"");
    
    // The input file/directory argument.
    addArgument(parser,
                ArgParseArgument(ArgParseArgument::INPUT_FILE, "IN"));

    // The output file argument.
    addOption(parser,
              ArgParseOption("o", "output-prefix", "output path prefix.",
                             ArgParseArgument::OUTPUT_FILE, "OUT"));
    addOption(parser,
              ArgParseOption("m",
                             "mapping-files",
                             "directory containing various mapping files "
                             "(gi2taxaID.map nodes_reduced.map "
                             "names_reduced.map taxaIDFeatures.map).",
                             ArgParseArgument::INPUT_FILE, "MAP-DIR"));
    
    addOption(parser,
              ArgParseOption("w",
                             "bin-width",
                             "Set the width of a single bin in neuclotides.",
                             ArgParseArgument::INTEGER, "INT"));
    addOption(parser,
              ArgParseOption("mr",
                             "min-reads",
                             "Minimum number of matching reads to consider "
                             "a reference present.",
                             ArgParseArgument::INTEGER, "INT"));
    addOption(parser,
              ArgParseOption("r", "rank",
                             "The taxonomic rank of identification", ArgParseOption::STRING));
    setValidValues(parser, "rank", options.rankList);
    setDefaultValue(parser, "rank", options.rank);
    
    
    setDefaultValue(parser, "bin-width", options.binWidth);
    setDefaultValue(parser, "min-reads", options.minReads);
 
    addOption(parser,
              ArgParseOption("c",
                             "cutoff",
                             "variation of the minnimum coverage "
                             "from mean (in standard deviation).",
                             ArgParseArgument::DOUBLE, "DOUBLE"));
    
    setMinValue(parser, "cutoff", "0");
    setMaxValue(parser, "cutoff", "3");
    setDefaultValue(parser, "cutoff", options.cutoff);
    
//    setDefaultValue(parser, "cov-cutoff", options.covCutoff);
//
//    addOption(parser,
//              ArgParseOption("ac",
//                            "abundance-cutoff",
//                             "minimum relative abundance of a taxon",
//                             ArgParseArgument::DOUBLE, "DOUBLE"));
//    
//    setDefaultValue(parser, "abundance-cutoff", options.abundanceCutoff);
 
    addOption(parser,
              ArgParseOption("d", "directory", "Input is a directory."));
    addOption(parser,
              ArgParseOption("v", "verbose", "Enable verbose output."));
    
    // Add Examples Section.
    addTextSection(parser, "Examples");

    addListItem(parser,
                "\\fBslimm\\fP \\fB-v\\fP \\fIexample.bam\\fP "
                "-o \\fIfeatures.tsv\\fP",
                "get the microbes present in \"example.bam\" alignment "
                "file and write the result to \"features.tsv\".");

    addListItem(parser, "\\fBslimm\\fP \\fB-d\\fP \\fIexample-dir\\fP "
                "-o \\fI~/results/\\fP",
                "get the microbes present in the SAM/BAM files located "
                "under \"example-dir\" and write the result to \""
                "in results folder under the home directory\".");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    
    if (res != ArgumentParser::PARSE_OK)
        return res;
    
    // Extract option values.
    if (isSet(parser, "mapping-files"))
        getOptionValue(options.mappingDir, parser, "mapping-files");

    if (isSet(parser, "bin-width"))
        getOptionValue(options.binWidth, parser, "bin-width");
    
    if (isSet(parser, "min-reads"))
        getOptionValue(options.minReads, parser, "min-reads");

    if (isSet(parser, "rank"))
        getOptionValue(options.rank, parser, "rank");
    
    if (isSet(parser, "cutoff"))
        getOptionValue(options.cutoff, parser, "cutoff");

    if (isSet(parser, "verbose"))
        getOptionValue(options.verbose, parser, "verbose");

    if (isSet(parser, "directory"))
        options.isDirectory = true;    

    getArgumentValue(options.inputPath, parser, 0);
    getOptionValue(options.outputPrefix, parser, "output-prefix");
    
    return ArgumentParser::PARSE_OK;
}
           
// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.
int main(int argc, char const ** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    
    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;
   
    // Prepare the
    StringList inputPaths;
    uint32_t numFiles = 1, totalRecCount = 0, fileCount = 0;
    if (options.isDirectory)
    {
        
        getFilesInDirectory(inputPaths, toCString(options.inputPath));
        numFiles = length(inputPaths);
        std::cout << numFiles << ": SAM/BAM Files found under the directory: "
        << options.inputPath <<"! \n" ;
    }
    else
    {
        if (is_file(toCString(options.inputPath)))
            appendValue(inputPaths, toCString(options.inputPath));
        else
        {
            std::cout << options.inputPath
            << " is not a file use -d option for a directory.\n";
            return 1;
        }
    }
    
    Iterator<StringList>::Type fileIt = begin(inputPaths);
    Iterator<StringList>::Type itEnd = end(inputPaths);
    
    std::stringstream ss;
    std::ofstream sam_extract_file;
    
    std::vector<std::string> uniqueReads;
    std::vector<std::string> uniqueReadsByTaxid;

    
    CharString giMap_path, features_path, nodes_path, names_path;
    giMap_path = features_path = nodes_path = names_path = options.mappingDir;
    
    append(giMap_path, "/gi2taxaID.map");
    append(features_path, "/taxaIDFeatures.map");
    append(nodes_path, "/nodes.dmp");
    append(names_path, "/names.dmp");

    Timer<> MainStopWatch;

    std::cout<<"Loding gi2taxaID mapping ... ";
    TIntIntMap gi2taxaID;
    
    gi2taxaID =  loadMapping<TIntIntMap, CharString>(giMap_path);
    std::cout<<"in " << MainStopWatch.lap() <<" secs [OK!]" <<std::endl;
    
    // ============get the taxaID2name mapping ================
    std::cout<<"Loding taxaID2name mapping ... ";
    TIntStrMap taxaID2name;
    
    taxaID2name =  loadMappingInt2String<TIntStrMap, CharString>(names_path);
    std::cout<<"in " << MainStopWatch.lap() <<" secs [OK!]" <<std::endl;
    
    // ============get the node mapping ================
    std::cout<<"Loding node mapping ... ";
    TNodes nodes;
    loadNodes<>(nodes, nodes_path);
    
    std::cout<<"in " << MainStopWatch.lap() <<" secs [OK!]" <<std::endl;
    
    // ============get the taxon properties ================
    std::cout<<"Loding additional taxon related features ...";
    std::unordered_map <int, TaxaProperty> taxonProperties;
    std::ifstream taxaIDMap(toCString(features_path));
    std::string line;
    while (std::getline(taxaIDMap, line))
    {
        std::vector<std::string> chunks = split(line, '\t');
        TaxaProperty tp;
        tp.name = chunks[1];
        tp.length = atof(toCString(chunks[2]));
        tp.GCContent = atof(toCString(chunks[3]));
        tp.neighborsLv1 = atoi(toCString(chunks[4]));
        tp.neighborsLv2 = atoi(toCString(chunks[5]));
        tp.neighborsLv3 = atoi(toCString(chunks[6]));
        
        taxonProperties[atoi(toCString(chunks[0]))]=tp;
    }
    taxaIDMap.close();
    std::cout<<"in " << MainStopWatch.lap() <<" secs [OK!]" <<std::endl;

    for(; fileIt != itEnd; goNext(fileIt))
    {
        Timer<> PerFileStopWatch;
        fileCount ++;
        CharString curr_file;
        
        std::cout<<"\n=============================\nReading "
        <<fileCount<<" of "<< numFiles<<" files ...\n"
                <<getFilename(toCString(value(fileIt)))<<"\n";
        // read the original file
        BamFileIn bamFile;
        if (!open(bamFile, toCString(value(fileIt))))
        {
            std::cerr << "Could not open " << toCString(value(fileIt)) << "!\n";
            return 1;
        }
        
        Slimm slimm;
        slimm.options = options;
        BamHeader header;
        readHeader(header, bamFile);
        
        StringSet<CharString> refNames = contigNames(context(bamFile));
        StringSet<uint32_t> refLengths;
        refLengths = contigLengths(context(bamFile));
        
        slimm.references.resize(length(refNames));
        uint32_t noOfRefs = length(refNames);
        std::cout<<"computing features of each reference genome ... ";

        for (uint32_t i=0; i<noOfRefs; ++i)
        {
            ReferenceContig current_ref;
            current_ref.refName = refNames[i];
            current_ref.length = refLengths[i];
            
            // Intialize coverages based on the length of a refSeq
            Coverage cov(current_ref.length, options.binWidth);
            current_ref.cov = cov;
            current_ref.uniqCov = cov;
            current_ref.uniqCov2 = cov;
            StringList chunks;
            strSplit(chunks, refNames[i], EqualsChar<'|'>());
            current_ref.taxaID = atoi(toCString(chunks[1]));
            
            slimm.matchedTaxa.push_back(current_ref.taxaID);
            slimm.references[i] = current_ref;
        }
        std::cout<<"in " << PerFileStopWatch.lap() <<" secs [OK!]" <<std::endl;
        
        
        std::cout<<"Analysing alignments, reads and references ...";
        analyzeAlignments(slimm, bamFile, options.binWidth);
        totalRecCount += slimm.noMatchedQueries;
        std::cout<<"in " << PerFileStopWatch.lap() <<" secs [OK!]" <<std::endl;
        
        float expCov = (slimm.avgQLength*slimm.noMatchedQueries) /
                        float(slimm.matchedRefsLen);

        std::cout<<"Expected Coverage = " << expCov <<std::endl;
        std::cout<<"Number of Ref with reads = " << slimm.noOfRefs <<std::endl;
        std::cout<<"Filtering unlikely sequences ..."  << std::endl ;
        filterAlignments(slimm);
        std::cout<<"in " << PerFileStopWatch.lap() <<" secs [OK!]" <<std::endl;

        std::cout   << length(slimm.validRefs)
                    <<" passed the threshould coverage."<< std::endl;
    
        std::cout<<"Writing features to a file ..." ;
        std::string outputFile = toCString(value(fileIt));
        if(length(options.outputPrefix) > 0)
        {
            outputFile = getDirectory(toCString(options.outputPrefix));
            outputFile.append("/");
            outputFile.append(getFilename(toCString(options.outputPrefix)));
            outputFile.append(getFilename(toCString(value(fileIt))));
        }
        std::string tsvFile = getTSVFileName(outputFile);
        writeToFile(tsvFile, slimm.references,
                    taxaID2name);

        std::cout<<"in " << PerFileStopWatch.lap() <<" secs [OK!]" <<std::endl;

        std::cout<<"Assigning reads to Least Common Ancestor (LCA) " ;
        getReadLCACount(slimm, nodes);
        
        tsvFile = getTSVFileName(outputFile, slimm.options.rank); 
        writeAbundance(slimm, nodes, taxaID2name, tsvFile) ;

        std::cout <<"in "
                  << PerFileStopWatch.lap() << " secs [OK!]" << std::endl;
        std::cout<<"File took " << PerFileStopWatch.elapsed()
        <<" secs to process." <<std::endl;
    }
    
    
    CharString output_directory = getDirectory(toCString(options.inputPath)) ;
    if(length(options.outputPrefix) > 0)
        output_directory = getDirectory(toCString(options.outputPrefix));

    std::cout << totalRecCount
              << " SAM/BAM alignment records are proccessed."<<std::endl;
    std::cout << "extracted features are written to: "
              << output_directory<<std::endl;
    std::cout << "Total tame elapsed: "
              << MainStopWatch.elapsed() <<" secs"<<std::endl;
    
    return 0;
}