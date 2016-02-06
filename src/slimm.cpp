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

#include <sys/stat.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <numeric>
#include <unordered_map>

#include "timer.h"
#include "load_mapping_files.h"


#include <seqan/basic.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

using namespace seqan;
typedef StringSet <Dna5String> SequenceList;
typedef StringSet <CharString> StringList;
typedef std::unordered_map <uint32_t, std::pair<uint32_t, std::string> > TNodes;


template <typename T1 = void, typename T2 = void>
struct Limits;

typedef StringSet <Dna5String>                              SequenceList;
typedef StringSet <CharString>                              StringList;

typedef __uint16                                            __intSizeGCount;
typedef __uint32                                            __intSizeBinWidth;
typedef __uint32                                            __intSizeQCount;
typedef __uint32                                            __intSizeGLength;
typedef __uint16                                            __intSizeQLength;
typedef __uint32                                            __intSizeMatchCount;
typedef StringSet<String<__intSizeMatchCount> >             TMatchSet;
typedef std::unordered_map <uint32_t, uint32_t>             TIntIntMap;
typedef std::unordered_map <uint32_t, float>                TIntFloatMap;
typedef std::unordered_map <uint32_t, std::string>          TIntStrMap;
typedef std::unordered_map <uint32_t, std::set<uint32_t> >  TIntSetMap;


// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------
using namespace seqan;
// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

inline uint32_t findLCATaxaID(std::set<uint32_t> const & taxaIDs,
                              TNodes const & nodes);

struct AppOptions
{
    __uint32            verbosity;  // 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    __uint32            cutoff = 3;  // 0 , 1 , 2 , 3.
    __intSizeBinWidth   binWidth = 100;
    bool                isDirectory;
    CharString          inputPath;
    CharString          outputPrefix;
    CharString          mappingDir;

    
    AppOptions() :
    verbosity(1),
    isDirectory(false),
    mappingDir("taxonomy/")
    {}
};


// ==========================================================================
// Classes
// ==========================================================================
class TaxaProperty
{
public:
    std::string name;
    unsigned int length;
    unsigned int neighborsLv1;
    unsigned int neighborsLv2;
    unsigned int neighborsLv3;
    float GCContent;
    TaxaProperty():
    name(""),
    length(0),
    neighborsLv1(0),
    neighborsLv2(0),
    neighborsLv3(0),
    GCContent(50.0)
    {}
};


class Coverage
{
public:
    CharString refName;
    uint32_t binWidth;
    bool hasNonZeroBins;
    uint32_t noOfBins;
    std::vector <uint32_t>  binsHeight;
    
    Coverage():
    binWidth(1000),
    hasNonZeroBins(false),
    noOfBins(0)
    {}
    
    Coverage(uint32_t totalLen, uint32_t width)
    {
        binWidth = width;
        noOfBins = totalLen/width + ((totalLen/width)*width < totalLen);
        hasNonZeroBins = false;
        std::vector<uint32_t> tmp (noOfBins, 0);
        binsHeight = tmp;
    }
    
};

class ReferenceContig
{
public:
    CharString refName;
    uint32_t length;
    uint32_t noOfReads;
    uint32_t noOfUniqReads;
    uint32_t noOfUniqReads2;
    uint32_t covDepth;
    Coverage cov;
    Coverage uniqCov;
    Coverage uniqCov2;
    uint32_t taxaID;
    float covPercent;
    float uniqCovPercent;
    float uniqCovPercent2;
    
    ReferenceContig():
    length(0),
    noOfReads(0),
    noOfUniqReads(0),
    noOfUniqReads2(0),
    covDepth(0),
    covPercent(0.0),
    uniqCovPercent(0.0),
    uniqCovPercent2(0.0)
    {}
};

enum class TaxnomicRank : uint16_t
{
    KINGDOM, PHYLUM , CLASS, ORDER, FAMILY, GENUS, SPECIES
};

class Read
{
public:
    std::vector<std::pair<__int32, uint32_t> > matchPoints;
    uint32_t sumRefLengths = 0;
    uint32_t len = 0;
    //checks if all the match points are in the same sequence
    bool isUniq();
    // checks if all the match points are in the same sequence
    // ignoring sequences that are not in refList
    bool isUniq(std::vector<uint32_t> const & taxaIDs,
                std::vector<uint32_t> const & valRefs);
    //checks if all the match points are in the same taxaID
    bool isUniq(std::vector<uint32_t> const & taxaIDs);
    void update(std::vector<uint32_t> const & refList,
                std::vector<ReferenceContig> const & references );
};

// ----------------------------------------------------------------------------
// Class Slimm
// ----------------------------------------------------------------------------
class Slimm
{
public:
    std::vector<ReferenceContig> references;
    std::vector<uint32_t> matchedRefs;
    std::vector<uint32_t> validRefs;
    std::map<CharString, Read> reads;
    __intSizeQLength avgQLength = 0;
    __intSizeGLength matchedRefsLen = 0;
    __intSizeMatchCount noOfRefs = 0;
    __intSizeQCount hitCount = 0;
    __intSizeQCount noMatchedQueries = 0;
    
    float covCutoff = 0.01;
    AppOptions options;
    
    TIntIntMap taxaID2ReadCount;
    TIntFloatMap taxaID2Abundance;
    TIntSetMap taxaID2Children;
};

//checks if all the match points are in the same sequence
bool Read::isUniq()
{
    size_t len = matchPoints.size();
    if (len == 0 || len == 1)
        return len;
    else
    {
        __int32 refID = matchPoints[0].first;
        for (size_t i=1; i < len; ++i)
            if(refID !=  matchPoints[i].first)
                return false;
    }
    return true;
}


// checks if all the match points are in the same sequence
// ignoring sequences that are not in refList
bool Read::isUniq(std::vector<uint32_t> const & taxaIDs,
                  std::vector<uint32_t> const & valRefs)
{
    size_t len = matchPoints.size();
    if (len == 0 || len == 1)
        return len;
    else
    {
        std::set<uint32_t> refTaxaIDs;
        for (size_t i=0; i < len; ++i)
        {
            uint32_t refID = taxaIDs[(matchPoints[i]).first];
            if(std::find(valRefs.begin(), valRefs.end(), refID) != valRefs.end())
                refTaxaIDs.insert(refID);
        }
        if (refTaxaIDs.size() > 1)
        {
            return false;
        }
    }
    return true;
}
//checks if all the match points are in the same taxaID
bool Read::isUniq(std::vector<uint32_t> const & taxaIDs)
{
    return isUniq(taxaIDs, taxaIDs);
}

void Read::update(std::vector<uint32_t> const & refList,
            std::vector<ReferenceContig> const & references )
{
//    if (length(matchPoints) == 0 || length(matchPoints) == 1)
    if (length(matchPoints) == 0)
        return;
    else
    {
        std::vector<std::pair<__int32, uint32_t> > newMatchPoints;
        std::vector<std::pair<__int32, uint32_t> >::iterator it;
        for (it=matchPoints.begin(); it!=matchPoints.end(); ++it)
        {
            if(!(std::find(refList.begin(), refList.end(),
                           it->first) == refList.end()))
            {
                newMatchPoints.push_back(*it);
            }
            else
            {
                sumRefLengths -= references[it->first].length;
            }
        }
        matchPoints = newMatchPoints;
    }
}



// ==========================================================================
// Functions
// ==========================================================================


bool is_file(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISREG(buf.st_mode);
}

bool is_dir(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISDIR(buf.st_mode);
}

// Calculates log2 of number.
float log_2(float n)
{
    // log(n)/log(2) is log2.
    return std::log(n)/std::log(2);
}

// returns a vector after spliting a string in two chunks
std::vector<std::string> &split(const std::string &s,
                                char delim,
                                std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

template <typename T>
std::string numberToString ( T Number )
{
    std::stringstream ss;
    ss<<Number;
    return ss.str();
}

template <typename T>
T stringToNumber ( const std::string &Text )
{
    std::stringstream ss(Text);
    T result;
    return ss >> result ? result : 0;
}

template <typename T>
T median(std::vector<T> &v)
{
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin()+n, v.end());
    if (n*2 == v.size())
        return v[n];
    else
        return (v[n] + v[n+1])/2.0;
}

template <typename T>
T mean(std::vector<T> &v) {
    T vSum = std::accumulate(v.begin(), v.end(), 0.0);
    return vSum/v.size();
}

template <typename T>
T variance(std::vector<T> &v, T m)
{
    T temp = 0;
    for(auto i : v)
    {
        temp += (i - m) * (i - m) ;
    }
    return temp / v.size();
}
    
template <typename T>
T variance(std::vector<T> &v)
{
    T m = mean(v);
    return variance(v, m);
}

template <typename T>
T stdDev(std::vector<T> &v, T m)
{
    return sqrt(variance(v, m));
}

template <typename T>
T stdDev(std::vector<T> &v)
{
    T m = mean(v);
    return sqrt(variance(v, m));
}


float calculateAlignmentScore(String<CigarElement<> > cigar,
                              int editDistance,
                              unsigned readLen)
{
    float score = 0.0;
    score += editDistance;
    typedef Iterator<String<CigarElement<> > >::Type TCigarIterator;
    
    for (TCigarIterator it = begin(cigar) ; it != end(cigar); goNext(it))
    {
        if (value(it).operation == 'D' || value(it).operation == 'I')
            score += float(value(it).count);
    }
    
    return score/readLen;
}

std::string getFilename (const std::string& str)
{
    std::size_t found = str.find_last_of("/\\");
    return str.substr(found+1);
}
std::string getDirectory (const std::string& str)
{
    std::size_t found = str.find_last_of("/\\");
    return str.substr(0,found);
}

std::string getTSVFileName (const std::string& fName)
{
    std::string result = fName;
    if(fName.find(".sam") == fName.find_last_of(".") ||
        fName.find(".bam")  == fName.find_last_of(".") )
        result.replace((result.find_last_of(".")), 4, "");
    result.append(".tsv");
    return result;
}
std::string getTSVFileName (const std::string& fName, const std::string& sfx)
{
    return (getTSVFileName(fName)).insert(fName.size()-4, sfx);
}


void getFilesInDirectory(StringList &inputPaths, std::string directory)
{
#ifdef WINDOWS
    HANDLE dir;
    WIN32_FIND_DATA file_data;
    
    if ((dir = FindFirstFile((directory + "/*").c_str(), &file_data)) == INVALID_HANDLE_VALUE)
        return; /* No files found */
    
    do {
        const std::string file_name = file_data.cFileName;
        const std::string full_file_name = directory + "/" + file_name;
        const bool is_directory = (file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;
        
        if (file_name[0] == '.')
            continue;
        
        if (is_directory)
            continue;
        
        if((full_file_name.find(".sam") == full_file_name.find_last_of(".") ||
            full_file_name.find(".bam")  == full_file_name.find_last_of(".") )
            appendValue(inputPaths, full_file_name);
        // out.push_back(full_file_name);
    } while (FindNextFile(dir, &file_data));
    
    FindClose(dir);
#else
    DIR *dir;
    struct dirent *ent;
    struct stat st;
    
    dir = opendir(directory.c_str());
    while ((ent = readdir(dir)) != NULL) {
        const std::string file_name = ent->d_name;
        const std::string full_file_name = directory + "/" + file_name;
        
        if (file_name[0] == '.')
            continue;
        
        if (stat(full_file_name.c_str(), &st) == -1)
            continue;
        
        const bool is_directory = (st.st_mode & S_IFDIR) != 0;
        
        if (is_directory)
            continue;
        
        
        if(full_file_name.find(".sam") == full_file_name.find_last_of(".") ||
           full_file_name.find(".bam")  == full_file_name.find_last_of(".") )
           appendValue(inputPaths, full_file_name);
        // out.push_back(full_file_name);
    }
    closedir(dir);
#endif
} // getFilesInDirectory

// ----------------------------------------------------------------------------
// Function setDateAndVersion()
// ----------------------------------------------------------------------------
   
void setDateAndVersion(ArgumentParser & parser)
{
    setCategory(parser, "Metagenomics");    
    #if defined(SEQAN_APP_VERSION) && defined(SEQAN_REVISION)
        setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    #endif
    #if defined(SEQAN_DATE)
        setDate(parser, SEQAN_DATE);
    #endif
}
   
// ----------------------------------------------------------------------------
// Function setDescription()
// ----------------------------------------------------------------------------
   
void setDescription(ArgumentParser & parser)
{
    addDescription(parser, "Species Level Identification of Microbes from Metagenomes");
    addDescription(parser, "See \\fIhttp://www.seqan.de/projects/slimm\\fP for more information.");
    addDescription(parser, "Investigates which microbial species are present from a BAM/SAM alignment file .");
    addDescription(parser, "(c) Copyright 2014-2017 by Temesgen H. Dadi.");
}
           
           
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
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "IN"));
    // The output file argument.
    addOption(parser, ArgParseOption("o", "output-prefix", "output path prefix.", 
                                            ArgParseArgument::OUTPUT_FILE, "OUT"));
    addOption(parser, ArgParseOption("m", "mapping-files", 
                                            "directory containing various mapping files "
                                            "(gi2taxaID.map nodes_reduced.map names_reduced.map taxaIDFeatures.map).",
                                            ArgParseArgument::INPUT_FILE, "MAP-DIR"));  
    addOption(parser, ArgParseOption("w", "bin-width", "Set the width of a single bin in neuclotides.",
                                            ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "bin-width", options.binWidth);
 
    addOption(parser, ArgParseOption("c", "cutoff", "variation of the minnimum coverage "
                                     "from mean in standard deviation.",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    setMinValue(parser, "cutoff", "0");
    setMaxValue(parser, "cutoff", "3");
    setDefaultValue(parser, "cutoff", options.cutoff);
    
//    setDefaultValue(parser, "cov-cutoff", options.covCutoff);
//
//    addOption(parser, ArgParseOption("ac", "abundance-cutoff", "minimum relative abundance of a taxon", 
//                                            ArgParseArgument::DOUBLE, "DOUBLE"));
//    setDefaultValue(parser, "abundance-cutoff", options.abundanceCutoff);
 
    addOption(parser, ArgParseOption("d", "directory", "Input is a directory."));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output.", 
                                            ArgParseArgument::INTEGER, "INT"));
    
    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBslimm\\fP \\fB-v\\fP \\fIexample.bam\\fP -o \\fIfeatures.tsv\\fP",
                "get the microbes present in \"example.bam\" alignment file and write the result to \"features.tsv\".");
    addListItem(parser, "\\fBslimm\\fP \\fB-d\\fP \\fIexample-dir\\fP -o \\fI~/results/\\fP",
                "get the microbes present in the SAM/BAM files located under \"example-dir\" and write the result to \""
                "in results folder under the home directory\".");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    
    if (res != ArgumentParser::PARSE_OK)
        return res;
    
    // Extract option values.
    if (isSet(parser, "mapping-files"))
        getOptionValue(options.mappingDir, parser, "mapping-files");
    if (isSet(parser, "bin-width"))
        getOptionValue(options.binWidth, parser, "bin-width");
    //    if (isSet(parser, "cov-cutoff"))
    //        getOptionValue(options.covCutoff, parser, "cov-cutoff");
    if (isSet(parser, "cutoff"))
        getOptionValue(options.cutoff, parser, "cutoff");
//    if (isSet(parser, "abundance-cutoff"))
//        getOptionValue(options.abundanceCutoff, parser, "abundance-cutoff");
    if (isSet(parser, "verbose"))
        getOptionValue(options.verbosity, parser, "verbose");
    if (isSet(parser, "directory"))
        options.isDirectory = true;    
    getArgumentValue(options.inputPath, parser, 0);
    getOptionValue(options.outputPrefix, parser, "output-prefix");
    
    return ArgumentParser::PARSE_OK;
}

           
// --------------------------------------------------------------------------
// Function getCovPercent()
// --------------------------------------------------------------------------
inline float getCovPercent(Coverage & cov)
{
    if(!cov.hasNonZeroBins)
        return 0.0;
    int nonZeroCount = 0;
    for (unsigned int i=0; i<cov.noOfBins; ++i )
        if (cov.binsHeight[i] > 0) {
            ++ nonZeroCount;
        }
    return float(nonZeroCount)/cov.noOfBins;
}
// --------------------------------------------------------------------------
// Function getCovDepth()
// --------------------------------------------------------------------------
inline float getCovDepth(Coverage & cov)
{
    if(!cov.hasNonZeroBins)
        return 0.0;
    
    std::vector <uint32_t>::iterator it;
    
    //copy the coverage height.
    std::vector <uint32_t>  bHeights = cov.binsHeight;
    std::vector <uint32_t>  bFrequentHeights;
    std::sort(bHeights.begin(), bHeights.end());
    
    uint32_t previousCount=1;
    uint32_t currentCount=1;
    
    for (unsigned int i=1; i<cov.noOfBins; ++i )
    {
        if (bHeights[i] == bHeights[i-1])
        {
            ++ currentCount;
        }
        else
        {
            if (currentCount==previousCount)
            {
                bFrequentHeights.push_back(bHeights[i-1]);
            }
            else if(currentCount > previousCount)
            {
                bFrequentHeights.empty();
                bFrequentHeights.push_back(bHeights[i-1]);
                previousCount = currentCount;
                
            }
            currentCount = 1;
        }
    }
    float sum = std::accumulate(bFrequentHeights.begin(), bFrequentHeights.end(), 0.0);
    return sum / length(bFrequentHeights);
}
uint32_t getLCA(std::vector<uint32_t> const & taxaIDs,
               std::vector<uint32_t> const & valTaxaIDs,
               TNodes const & nodes)
{
   //consider only those under validTaxaIDs
   std::set<uint32_t> parents;
   for (auto tID : taxaIDs)
   {
       if (std::find(valTaxaIDs.begin(), valTaxaIDs.end(), tID) != valTaxaIDs.end())
           parents.insert(tID);
   }
   while (parents.size() > 1)
   {
       std::set<uint32_t> newParents;
       uint32_t current_count = parents.size();
       for (std::set<uint32_t>::iterator it = parents.begin();
            it!= parents.end(); ++it)
       {
           if (nodes.find(*it) != nodes.end())
           {
               if(parents.find((nodes.at(*it).first)) != parents.end())
               {
                   --current_count;
                   if(current_count == 1)
                       return (nodes.at(*it)).first;
                   else
                       continue;
               }
               else
                   newParents.insert(nodes.at(*it).first);
           }
       }
       
       parents = newParents;
   }
   return *(parents.begin());
}
           
uint32_t getLCA(std::vector<uint32_t> const & taxaIDs,
               TNodes const & nodes)
{
    return getLCA(taxaIDs, taxaIDs, nodes);
}
           
inline void analyzeAlignments(Slimm & slimm,
                              BamFileIn & bamFile,
                              __intSizeBinWidth const & binWidth)
{
    BamAlignmentRecord record;
    while (!atEnd(bamFile))
    {
        readRecord(record, bamFile);
        __intSizeQLength readLen = record._l_qseq;
//        if (hasFlagUnmapped(record) || record.rID == BamAlignmentRecord::INVALID_REFID)
        if (record.rID == BamAlignmentRecord::INVALID_REFID)
            continue;  // Skip these records.
        
        uint32_t relativeBinNo = (record.beginPos + (readLen/2))/binWidth;
        slimm.references[record.rID].cov.hasNonZeroBins = true;
        
        // maintain read properties under slimm.reads
        CharString readName = record.qName;
        if(hasFlagFirst(record))
            append(readName, ".1");
        else if(hasFlagLast(record))
            append(readName, ".2");
        if (slimm.reads.count(readName) == 1)
        {
            slimm.reads[readName].matchPoints.push_back(std::pair<__int32, uint32_t>(record.rID, relativeBinNo));
        }
        else
        {
            Read query;
            query.matchPoints.push_back(std::pair<__int32, uint32_t>(record.rID, relativeBinNo));
            query.len = readLen;
            slimm.reads.insert(std::pair<CharString, Read>(readName, query));
        }
        ++slimm.hitCount;
    }
    slimm.noMatchedQueries = length(slimm.reads);
    unsigned totalUniqueReads = 0;
    __intSizeGLength conctQLength = 0;
    std::map<CharString, Read>::iterator it;
    for (it= slimm.reads.begin(); it != slimm.reads.end(); ++it)
    {
        conctQLength += it->second.len;
        if(it->second.isUniq(slimm.matchedRefs))
        {
            __int32 rID = it->second.matchPoints[0].first;
            uint32_t binNo = it->second.matchPoints[0].second;
            ++slimm.references[rID].noOfUniqReads;
            slimm.references[rID].uniqCov.hasNonZeroBins = true;
            ++slimm.references[rID].uniqCov.binsHeight[binNo];
            ++totalUniqueReads;
            ++slimm.references[rID].noOfReads;
            it->second.sumRefLengths += slimm.references[record.rID].length;
        }
        else
        {
            std::set<__int32> matchedRefs;
            std::for_each(it->second.matchPoints.begin(),
                          it->second.matchPoints.end(),
                          [&](std::pair<__int32, uint32_t> &mp)
            {
                size_t oldSize = matchedRefs.size();
                matchedRefs.insert(mp.first);
                size_t newSize = matchedRefs.size();
                if (newSize > oldSize)
                {
                    ++slimm.references[mp.first].noOfReads;
                    slimm.references[mp.first].cov.hasNonZeroBins = true;
                    ++slimm.references[mp.first].cov.binsHeight[mp.second];
                    it->second.sumRefLengths += slimm.references[mp.first].length;
                }
            });
            
        }
    }
    std::vector<float> covValues;
    slimm.avgQLength = conctQLength/slimm.noMatchedQueries;
    for (uint32_t i=0; i<length(slimm.references); ++i)
    {
        if (slimm.references[i].noOfReads > 0)
        {
            ++slimm.noOfRefs;
            slimm.matchedRefsLen += slimm.references[i].length;
            slimm.references[i].covDepth = getCovDepth(slimm.references[i].cov);
            slimm.references[i].covPercent = getCovPercent(slimm.references[i].cov);
            if(slimm.references[i].covPercent > 0.0)
                covValues.push_back(log(slimm.references[i].covPercent));
            else
                continue;
            slimm.references[i].uniqCovPercent = getCovPercent(slimm.references[i].uniqCov);
        }
    }
    float m = mean(covValues);
    float sd = stdDev(covValues, m);
    slimm.covCutoff = exp(m - slimm.options.cutoff *sd);
    std::cout<<"Mean = " << m <<" SD = " << sd <<" Cutoff = " << slimm.covCutoff <<std::endl;
    std::cout<< slimm.hitCount <<" recoreds processed."<<std::endl;
    std::cout<< "\t" << slimm.noMatchedQueries << " matching reads" <<std::endl;
    std::cout<< "\t" << totalUniqueReads <<" uniquily matching reads"<<std::endl;
   
}
           
inline void filterAlignments(Slimm & slimm)
{
    uint32_t noOfRefs = length(slimm.references);
    for (uint32_t i=0; i<noOfRefs; ++i)
    {
        if (slimm.references[i].noOfReads == 0)
            continue;
        if (slimm.references[i].covPercent > slimm.covCutoff &&
            slimm.references[i].noOfReads > 100)
            slimm.validRefs.push_back(i);
        else
            slimm.references[i]=ReferenceContig();
        
    }

    std::map<CharString, Read>::iterator it;
    uint32_t totalUniqueReads = 0;
    for (it= slimm.reads.begin(); it != slimm.reads.end(); ++it)
    {
        it->second.update(slimm.validRefs, slimm.references);
        if(it->second.isUniq(slimm.matchedRefs, slimm.validRefs))
        {
            __int32 rID = it->second.matchPoints[0].first;
            uint32_t binNo = it->second.matchPoints[0].second;
            slimm.references[rID].uniqCov2.hasNonZeroBins = true;
            ++slimm.references[rID].noOfUniqReads2;
            ++slimm.references[rID].uniqCov2.binsHeight[binNo];
            ++totalUniqueReads;
        }
    }
    std::cout << "Total number of uniquily matching reads "
    "(after recomputition) = " << totalUniqueReads <<std::endl;
}
           
           
           
           
           
inline void writeToFile(std::string & filePath,
                        std::vector<ReferenceContig> & refList,
//                        std::unordered_map <int, TaxaProperty> & taxonProperties,
                        TIntStrMap const & taxaID2name)
{
    std::ofstream features_file;
    features_file.open(filePath);
    
    features_file <<    "No.\t"
                        "Candidate Name\t"
                        "Taxid\t"
                        "Gen. Length\t"
                        "Num. Reads\t"
                        "Unique Reads\t"
                        "Unique Reads2\t"
                        "Coverage Depth\t"
                        "Mapping Error\t"
                        "Unique read cov_percentage\t"
                        "Unique read cov_percentage2\n";
    
    uint32_t current = 0;
    uint32_t noOfRefs = length(refList);
    for (uint32_t i=0; i < noOfRefs; ++i)
    {
        current ++;
        refList[i].uniqCovPercent2 = getCovPercent(refList[i].uniqCov2);
        ReferenceContig current_ref = refList[i];
        CharString candidateName = current_ref.refName;
        TIntStrMap::const_iterator it = taxaID2name.find(current_ref.taxaID);
        if (it != taxaID2name.end())
            candidateName = (taxaID2name.at(current_ref.taxaID));
        
        features_file   << current << "\t"
                        << candidateName << "\t"
                        << current_ref.taxaID << "\t"
                        << current_ref.length << "\t"
                        << current_ref.noOfReads << "\t"
                        << current_ref.noOfUniqReads << "\t"
                        << current_ref.noOfUniqReads2 << "\t"
                        << current_ref.covDepth << "\t"
                        << current_ref.uniqCovPercent << "\t"
                        << current_ref.uniqCovPercent2 << "\n";
    }
    features_file.close();
}

inline void getReadLCACount(Slimm & slimm,
                           TNodes const & nodes)
{
   std::map<CharString, Read>::iterator it;
   for (it= slimm.reads.begin(); it != slimm.reads.end(); ++it)
   {
       if (it->second.sumRefLengths > 0)
       {
           std::vector<uint32_t> tIDs;
           std::vector<uint32_t> rIDs;
           std::vector<std::pair<__int32, uint32_t> >::iterator it2;
           for (it2 = (it->second.matchPoints).begin();
                it2 != (it->second.matchPoints).end(); ++it2)
           {
               tIDs.push_back(slimm.matchedRefs[it2->first]);
               rIDs.push_back(it2->first);
           }
           uint32_t lcaTaxaID = getLCA(tIDs, nodes);
           // If taxaID already exists
           if (slimm.taxaID2ReadCount.count(lcaTaxaID) == 1)
               ++slimm.taxaID2ReadCount[lcaTaxaID];
           else   // first time for taxaID
               slimm.taxaID2ReadCount[lcaTaxaID] = 1;
           //add the contributing children
           for(auto i : rIDs)
               slimm.taxaID2Children[lcaTaxaID].insert(i);
           //add the read to all ancestors of the LCA
           while (nodes.count(lcaTaxaID) == 1 && lcaTaxaID != 0)
           {
               lcaTaxaID = (nodes.at(lcaTaxaID)).first;
               if (slimm.taxaID2ReadCount.count(lcaTaxaID) == 1)
                   ++slimm.taxaID2ReadCount[lcaTaxaID];
               else
                   slimm.taxaID2ReadCount[lcaTaxaID] = 1;
               //add the contributing children
               for(auto i : rIDs)
                   slimm.taxaID2Children[lcaTaxaID].insert(i);
           }
       }
   }
}

//inline void writeAbundance2(std::string const & filePath,
//                           TIntIntMap const & taxaIDReadCount,
//                           TNodes & nodes, TIntStrMap const & taxaID2name,
//                           std::string rank,
//                           uint32_t const & totalReads,
//                           float const & cutoff)
//{
//    std::ofstream abundunceFile;
//    abundunceFile.open(filePath);
//    abundunceFile<<"No.\tName\tTaxid\tNoOfReads\tRelativeAbundance\n";
//    
//    // calculate the total number of reads matching uniquily at given rank.
//    uint32_t noReadsAtRank = 0;
//    for (auto tID : taxaIDReadCount) {
//        if (rank == nodes[tID.first].second)
//            noReadsAtRank +=  tID.second ;
//    }
//    uint32_t unknownReads = totalReads-noReadsAtRank;
//    uint32_t count = 1;
//    for (auto tID : taxaIDReadCount) {
//        
//        if (rank == nodes[tID.first].second)
//        {
//            float abundance = float(tID.second)/totalReads;
//            // If the abundance is lower than a threshold do not report it
//            // Put the reads under the unkown
//            if (abundance < cutoff)
//            {
//                unknownReads += tID.second;
//                continue;
//            }
//            TIntStrMap::const_iterator it = taxaID2name.find (tID.first);
//            CharString candidateName = "Organism name not found";
//            if (it != taxaID2name.end())
//                candidateName = (taxaID2name.at(tID.first));
//            abundunceFile << count << "\t" << candidateName << "\t" << tID.first << "\t"
//                        << tID.second << "\t"  << abundance << "\n";
//            ++count;
//        }
//    }
//    // add the remaining  matching reads to unknowns.
//    abundunceFile   << count << "\tunknown_"<<rank<< "(multiple)" << "\t0\t"
//                    << unknownReads << "\t" << float(unknownReads)/totalReads << "\n";
//    abundunceFile.close();
//    
//
//}

//inline void writeAbundance(std::string const & filePath,
//                           TIntIntMap const & taxaIDReadCount,
//                           TIntFloatMap const & taxaIDAbundance,
//                           TNodes & nodes, TIntStrMap const & taxaID2name,
//                           std::string rank,
//                           uint32_t const & totalReads,
//                           float const & cutoff)
//{
//    std::ofstream abundunceFile;
//    abundunceFile.open(filePath);
//    
//    
//    abundunceFile<<"No.\tName\tTaxid\tNoOfReads\tRelativeAbundance\n";
//    
//    // calculate the total number of reads matching uniquily at that species level.
//    uint32_t noReadsAtRank = 0;
//    float abundancesAtRank = 0;
//    for (auto tID : taxaIDReadCount) {
//        if (rank == nodes[tID.first].second)
//        {
//            noReadsAtRank +=  tID.second ;
//            abundancesAtRank += taxaIDAbundance.at(tID.first);
//        }
//    }
//    uint32_t unknownReads = totalReads-noReadsAtRank;
//    float unknownAbundance = 0.0;
//    uint32_t count = 1;
//    for (auto tID : taxaIDReadCount) {
//        
//        if (rank == nodes[tID.first].second)
//        {
//            float abundance = taxaIDAbundance.at(tID.first)/abundancesAtRank;
//            // If the abundance is lower than a threshold do not report it
//            // Put the reads under the unkown
//            if (abundance < cutoff)
//            {
//                unknownReads += tID.second;
//                unknownAbundance += abundance;
//                continue;
//            }
//            std::unordered_map <uint32_t, std::string>::const_iterator it = taxaID2name.find (tID.first);
//            seqan::CharString candidateName = "Organism name not found";
//            if (it != taxaID2name.end())
//                candidateName = (taxaID2name.at(tID.first));
//            abundunceFile << count << "\t" << candidateName << "\t" << tID.first << "\t"
//            << tID.second << "\t"  << abundance << "\n";
//            ++count;
//        }
//    }
//    // add the remaining  matching reads to unknowns.
//    abundunceFile   << count << "\tunknown_"<<rank<< "(multiple)" << "\t0\t"
//    << unknownReads << "\t" << unknownAbundance << "\n";
//    abundunceFile.close();
//}

//inline void writeAbundance(Slimm const & slimm,
//                              TNodes & nodes, TIntStrMap const & taxaID2name,
//                              std::string const & filePath,
//                              std::string rank)
//{
//    std::ofstream abundunceFile;
//    abundunceFile.open(filePath);
//    
//    
//    abundunceFile<<"No.\tName\tTaxid\tNoOfReads\tCoverage\tRelativeAbundance\n";
//    
//    // calculate the total number of reads matching uniquily at that species level.
//    uint32_t noReadsAtRank = 0;
//    for (auto tID : slimm.taxaID2ReadCount) {
//        if (rank == nodes[tID.first].second)
//        {
//            noReadsAtRank +=  tID.second ;
//        }
//    }
//    uint32_t unknownReads = slimm.noMatchedQueries-noReadsAtRank;
//    float unknownAbundance = 0.0;
//    uint32_t count = 1;
//    uint32_t faild_count = 0;
//    for (auto tID : slimm.taxaID2ReadCount) {
//        if (rank == nodes[tID.first].second)
//        {
//            uint32_t contributersLength = 0;
//            std::set<uint32_t>::iterator it;
//            for (it=slimm.taxaID2Children.at(tID.first).begin();
//                 it!=slimm.taxaID2Children.at(tID.first).end(); ++it)
//                contributersLength += slimm.references[*it].length;
//            float factor = slimm.matchedRefsLen / contributersLength;
//            float coverage = float(tID.second * slimm.avgQLength)/contributersLength;
//            float relAbundance = float(tID.second * factor)/noReadsAtRank;
//            // If the abundance is lower than a threshold do not report it
//            // Put the reads under the unkown
//            if (relAbundance < slimm.abundanceCutoff)
//            {
//                unknownReads += tID.second;
//                unknownAbundance += relAbundance;
//                ++faild_count;
//                continue;
//            }
//            std::unordered_map <uint32_t, std::string>::const_iterator it2 =
//            taxaID2name.find (tID.first);
//            seqan::CharString candidateName = "Organism name not found";
//            if (it2 != taxaID2name.end())
//                candidateName = (taxaID2name.at(tID.first));
//            abundunceFile << count << "\t" << candidateName << "\t" << tID.first << "\t"
//            << tID.second << "\t"  << coverage << "\t"  << relAbundance << "\n";
//            ++count;
//        }
//    }
//    
//    // add the remaining  matching reads to unknowns.
//    abundunceFile   << count << "\tunknown_"<<rank<< "(multiple)" << "\t0\t"
//    << unknownReads << "\t0\t" << unknownAbundance << "\n";
//    abundunceFile.close();
//    std::cout<< faild_count <<" bellow cutoff ("<< slimm.abundanceCutoff <<") ...";
//}
   
           
inline void writeAbundance(Slimm const & slimm,
                           TNodes & nodes, TIntStrMap const & taxaID2name,
                           std::string const & filePath,
                           std::string rank)
{
    std::ofstream abundunceFile;
    abundunceFile.open(filePath);


    abundunceFile<<"No.\tName\tTaxid\tNoOfReads\tCoverage\tRelativeAbundance\n";

    // calculate the total number of reads matching uniquily at that species level.
    uint32_t noReadsAtRank = 0;
    for (auto tID : slimm.taxaID2ReadCount) {
        if (rank == nodes[tID.first].second)
        {
            noReadsAtRank +=  tID.second ;
        }
    }
    uint32_t unknownReads = slimm.noMatchedQueries-noReadsAtRank;
    float unknownAbundance = 0.0;
    uint32_t count = 1;
    uint32_t faild_count = 0;
    TIntFloatMap cladeCov;
    float totalCov = 0.0;
    std::vector<float> covValues;
    for (auto tID : slimm.taxaID2ReadCount) {
        if (rank == nodes[tID.first].second)
        {
            uint32_t contributersLength = 0;
            std::set<uint32_t>::iterator it;
            for (it=slimm.taxaID2Children.at(tID.first).begin();
                 it!=slimm.taxaID2Children.at(tID.first).end(); ++it)
                contributersLength += slimm.references[*it].length;
            cladeCov[tID.first] = float(tID.second * slimm.avgQLength)/contributersLength;
            if(cladeCov[tID.first] > 0.0)
                covValues.push_back(log(cladeCov[tID.first]));
            totalCov += cladeCov[tID.first];
        }
    }

    float m = mean(covValues);
    float sd = stdDev(covValues, m);
    float cutoff = exp(m - slimm.options.cutoff*sd);
    std::cout<<"Mean = " << m <<" SD = " << sd <<" Cutoff = " << cutoff <<std::endl;
    
    for (auto tID : cladeCov) {
        float relAbundance = tID.second/totalCov;
        // If the abundance is lower than a threshold do not report it
        // Put the reads under the unkown
        if (relAbundance < 0.0)
        {
            unknownReads += totalCov;
            unknownAbundance += relAbundance;
            ++faild_count;
            continue;
        }
        std::unordered_map <uint32_t, std::string>::const_iterator it2 =
        taxaID2name.find (tID.first);
        seqan::CharString candidateName = "Organism name not found";
        if (it2 != taxaID2name.end())
            candidateName = (taxaID2name.at(tID.first));
        abundunceFile << count << "\t" << candidateName << "\t" << tID.first << "\t"
        << slimm.taxaID2ReadCount.at(tID.first) << "\t"  << tID.second << "\t"  << relAbundance << "\n";
        ++count;
    }
    
    // add the remaining  matching reads to unknowns.
    abundunceFile   << count << "\tunknown_"<<rank<< "(multiple)" << "\t0\t"
    << unknownReads << "\t0\t" << unknownAbundance << "\n";
    abundunceFile.close();
    std::cout<< faild_count <<" bellow cutoff ("<< cutoff <<") ...";
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
            
            slimm.matchedRefs.push_back(current_ref.taxaID);
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

        std::cout<< length(slimm.validRefs) <<" passed the threshould coverage."<< std::endl;
    
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
//                    taxonProperties,
                    taxaID2name);

        std::cout<<"in " << PerFileStopWatch.lap() <<" secs [OK!]" <<std::endl;

        std::cout<<"assigning reads to Least Common Ancestor (LCA) " ;
        getReadLCACount(slimm, nodes);
        
        tsvFile = getTSVFileName(outputFile, "_sp_reported");
        writeAbundance(slimm, nodes, taxaID2name, tsvFile, "species") ;

        std::cout <<"in " << PerFileStopWatch.lap() << " secs [OK!]" << std::endl;
        std::cout<<"File took " << PerFileStopWatch.elapsed()
        <<" secs to process." <<std::endl;
    }
    
    
    CharString output_directory = getDirectory(toCString(options.inputPath)) ;
    if(length(options.outputPrefix) > 0)
        output_directory = getDirectory(toCString(options.outputPrefix));
    std::cout<<totalRecCount<<" SAM/BAM alignment records are proccessed."<<std::endl;
    std::cout<<"extracted features are written to: " << output_directory<<std::endl;
    std::cout<<"Total tame elapsed: " << MainStopWatch.elapsed() <<" secs"<<std::endl;
    
    return 0;
}
