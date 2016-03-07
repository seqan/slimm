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

#ifndef slimm_h
#define slimm_h

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
typedef std::pair<__int32, uint32_t>                        TMatchPoint;
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
    __uint32            cutoff = 3;  // 0 , 1 , 2 , 3.
    __intSizeBinWidth   binWidth = 100;
    __intSizeBinWidth   minReads = 100;
    bool                verbose;
    bool                isDirectory;
    CharString          inputPath;
    CharString          outputPrefix;
    CharString          mappingDir;
    
    
    AppOptions() :
    verbose(false),
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
    std::string         name;
    unsigned int        length;
    unsigned int        neighborsLv1;
    unsigned int        neighborsLv2;
    unsigned int        neighborsLv3;
    float               GCContent;
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
    CharString                  refName;
    uint32_t                    binWidth;
    uint32_t                    noOfBins;
    uint32_t                    noOfNonZeroBins;
    std::vector <uint32_t>      binsHeight;
    
    Coverage():
    binWidth(1000),
    noOfBins(0),
    noOfNonZeroBins(0)
    {}
    
    Coverage(uint32_t totalLen, uint32_t width)
    {
        binWidth = width;
        noOfBins = totalLen/width + ((totalLen/width)*width < totalLen);
        noOfNonZeroBins = 0;
        std::vector<uint32_t> tmp (noOfBins, 0);
        binsHeight = tmp;
    }
    
};

class ReferenceContig
{
public:
    CharString          refName;
    uint32_t            length;
    uint32_t            noOfReads;
    uint32_t            noOfUniqReads;
    uint32_t            noOfUniqReads2;
    Coverage            cov;
    Coverage            uniqCov;
    Coverage            uniqCov2;
    uint32_t            taxaID;
    
    ReferenceContig():
    length(0),
    noOfReads(0),
    noOfUniqReads(0)
    {}
    
    float               uniqCovPercent();
    float               uniqCovPercent2();
    uint32_t            covDepth();
    float               covPercent();
};


class TargetRef
{
public:
    __int32                     rID;
    std::vector<uint32_t>      positions;

    //constructer takes a ref id and a position for the first time
    TargetRef(__int32 ref, uint32_t pos)
    {
        rID = ref;
        positions.push_back(pos);
    }
};


enum class TaxnomicRank : uint16_t
{
    KINGDOM, PHYLUM , CLASS, ORDER, FAMILY, GENUS, SPECIES
};

class Read
{
public:
    std::vector<TargetRef>          targets;
    uint32_t                        sumRefLengths = 0;
    uint32_t                        len = 0;
    
    
    //checks if all the match points are in the same sequence
    bool isUniq();
    
    // checks if all the match points are in the same sequence
    // ignoring sequences that are not in refList
    bool isUniq(std::vector<uint32_t> const & taxaIDs,
                std::set<uint32_t> const & valRefs);
    
    //checks if all the match points are in the same taxaID
    bool isUniq(std::vector<uint32_t> const & taxaIDs);
    
    void update(std::set<uint32_t> const & valRefs,
                std::vector<ReferenceContig> const & references );
    
    void addTarget(int32_t rID, uint32_t binNo);

};

// ----------------------------------------------------------------------------
// Class Slimm
// ----------------------------------------------------------------------------
class Slimm
{
public:
    std::vector<ReferenceContig>    references;
    std::vector<uint32_t>           matchedTaxa;
    std::set<uint32_t>              validRefs;
    std::unordered_map<std::string, Read>      reads;
    __intSizeQLength                avgQLength = 0;
    __intSizeGLength                matchedRefsLen = 0;
    __intSizeMatchCount             noOfRefs = 0;
    __intSizeQCount                 hitCount = 0;
    __intSizeQCount                 noMatchedQueries = 0;
    
    AppOptions                      options;
    
    TIntIntMap                      taxaID2ReadCount;
    TIntFloatMap                    taxaID2Abundance;
    TIntSetMap                      taxaID2Children;
    float                           covCutoff = 0.01;
};

//checks if all the match points are in the same sequence
bool Read::isUniq()
{
    size_t len = targets.size();
    if (len == 0 || len == 1)
        return len;
    return false;
}


// checks if all the match points are in the same sequence
// ignoring sequences that are not in refList
bool Read::isUniq(std::vector<uint32_t> const & taxaIDs,
                  std::set<uint32_t> const & valRefs)
{
    size_t len = targets.size();
    if (len == 0 || len == 1)
        return len;
    else
    {
        std::set<uint32_t> refTaxaIDs;
        for (size_t i=0; i < len; ++i)
        {
            uint32_t refID = taxaIDs[(targets[i]).rID];
            if(valRefs.find(refID) != valRefs.end())
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
//bool Read::isUniq(std::vector<uint32_t> const & taxaIDs)
//{
//    std::set<uint32_t> s(taxaIDs.begin(), taxaIDs.end());
//    return isUniq(taxaIDs, s);
//}

void Read::update(std::set<uint32_t> const & valRefs,
                  std::vector<ReferenceContig> const & references )
{
    size_t len = targets.size();
    if (len == 0 || len == 1)
        return;
    else
    {
        std::vector<TargetRef> newTargets;
        for (size_t i=0; i < len; ++i)
        {
            if(valRefs.find((targets[i]).rID) != valRefs.end())
                newTargets.push_back(targets[i]);
            else
                sumRefLengths -= references[(targets[i]).rID].length;
        }
        targets = newTargets;
    }
}

void Read::addTarget(int32_t rID, uint32_t binNo)
{
    size_t len = targets.size();
    if (len == 0 )
    {
        targets.push_back(TargetRef(rID, binNo));
        return;
    }
    else
    {
        for (size_t i=0; i < len; ++i)
        {
            if((targets[i]).rID == rID)
            {
                targets[i].positions.push_back(binNo);
                return;
            }
        }
        targets.push_back(TargetRef(rID, binNo));
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

// returns a vector after spliting a string into two chunks
std::vector<std::string> &split(const std::string &s,
                                char delim,
                                std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}
std::vector<std::string> split(const std::string &s, char delim)
{
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
    addDescription(parser,
                   "Species Level Identification of Microbes from Metagenomes");
    addDescription(parser,
                   "See \\fIhttp://www.seqan.de/projects/slimm\\fP "
                   "for more information.");
    addDescription(parser,
                   "Investigates which microbial species are present from "
                   "a BAM/SAM alignment file .");
    addDescription(parser, "(c) Copyright 2014-2017 by Temesgen H. Dadi.");
}

// --------------------------------------------------------------------------
// Function getCovPercent()
// --------------------------------------------------------------------------
inline float getCovPercent(Coverage cov)
{
    return float(cov.noOfNonZeroBins)/cov.noOfBins;
}
float ReferenceContig::covPercent()
{
    return getCovPercent(cov);
}
float ReferenceContig::uniqCovPercent()
{
    return getCovPercent(uniqCov);
}
float ReferenceContig::uniqCovPercent2()
{
    return getCovPercent(uniqCov2);
}

// --------------------------------------------------------------------------
// Function getCovDepth()
// --------------------------------------------------------------------------
uint32_t ReferenceContig::covDepth()
{
    if(cov.noOfNonZeroBins == 0)
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
    float sum = std::accumulate(bFrequentHeights.begin(),
                                bFrequentHeights.end(), 0.0);
    return sum / bFrequentHeights.size();
}
uint32_t getLCA(std::set<uint32_t> const & taxaIDs,
                std::set<uint32_t> const & valTaxaIDs,
                TNodes const & nodes)
{
    //consider only those under validTaxaIDs
    std::set<uint32_t> parents;
    for (auto tID : taxaIDs)
    {
        if (valTaxaIDs.find(tID) != valTaxaIDs.end())
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

uint32_t getLCA(std::set<uint32_t> const & taxaIDs,
                TNodes const & nodes)
{
    return getLCA(taxaIDs, taxaIDs, nodes);
}

uint32_t getLCA(std::vector<uint32_t> const & taxaIDs,
                TNodes const & nodes)
{
    std::set<uint32_t> s(taxaIDs.begin(), taxaIDs.end());
    return getLCA(s, s, nodes);
}


inline void analyzeAlignments(Slimm & slimm,
                              BamFileIn & bamFile,
                              __intSizeBinWidth const & binWidth)
{

    BamAlignmentRecord record;
    while (!atEnd(bamFile))
    {
        readRecord(record, bamFile);
        if (hasFlagUnmapped(record) || record.rID == BamAlignmentRecord::INVALID_REFID)
            continue;  // Skip these records.
        
        uint32_t queryLen = length(record.seq);
        uint32_t relativeBinNo = (record.beginPos + (queryLen/2))/binWidth;
        ++slimm.references[record.rID].cov.noOfNonZeroBins;
        
        // maintain read properties under slimm.reads
        std::string readName = toCString(record.qName);
        if(hasFlagFirst(record))
            append(readName, ".1");
        else if(hasFlagLast(record))
            append(readName, ".2");
//      if there is no read with readName this will create one.
        slimm.reads[readName].addTarget(record.rID, relativeBinNo);
        slimm.reads[readName].len = queryLen;
        ++slimm.hitCount;
    }
    unsigned totalUniqueReads =0;
    __intSizeGLength conctQLength = 0;

    std::set<uint32_t> s(slimm.matchedTaxa.begin(), slimm.matchedTaxa.end());
    
    for (auto it= slimm.reads.begin(); it != slimm.reads.end(); ++it)
    {
        conctQLength += it->second.len;
        if(it->second.isUniq(slimm.matchedTaxa, s))
        {
            __int32 rID = it->second.targets[0].rID;
            uint32_t binNo = it->second.targets[0].positions[0];
            ++slimm.references[rID].noOfUniqReads;
            ++slimm.references[rID].uniqCov.noOfNonZeroBins;
            ++slimm.references[rID].uniqCov.binsHeight[binNo];
            ++totalUniqueReads;
            ++slimm.references[rID].noOfReads;
            it->second.sumRefLengths += slimm.references[record.rID].length;
        }
        else
        {
            size_t len = it->second.targets.size();
            for (size_t i=0; i < len; ++i)
            {

                __int32 refID = (it->second.targets[i]).rID;
                uint32_t p = (it->second.targets[i]).positions[0];
                
                //only the first match will be counted
                
                ++slimm.references[refID].noOfReads;
                ++slimm.references[refID].cov.noOfNonZeroBins;
                ++slimm.references[refID].cov.binsHeight[p];
                it->second.sumRefLengths += slimm.references[refID].length;
            }
        }
    }
    slimm.noMatchedQueries = slimm.reads.size();

    std::vector<float> covValues;
    slimm.avgQLength = conctQLength/slimm.noMatchedQueries;
    for (uint32_t i=0; i<length(slimm.references); ++i)
    {
        if (slimm.references[i].noOfReads > 0)
        {
            ++slimm.noOfRefs;
            slimm.matchedRefsLen += slimm.references[i].length;
            if(slimm.references[i].covPercent() > 0.0)
                covValues.push_back(log(slimm.references[i].covPercent()));
            else
                continue;
        }
    }
    float m = mean(covValues);
    float sd = stdDev(covValues, m);
    slimm.covCutoff = exp(m - slimm.options.cutoff *sd);
    std::cout   << std::endl
                << "Mean = " << m
                <<" SD = " << sd
                <<" Cutoff = " << slimm.covCutoff <<std::endl;
    
    std::cout << slimm.hitCount << " recoreds processed." << std::endl;
    std::cout << "\t" << slimm.noMatchedQueries << " matching reads" << std::endl;
    std::cout << "\t" << totalUniqueReads << " uniquily matching reads"<< std::endl;
    
}

inline void filterAlignments(Slimm & slimm)
{
    uint32_t noOfRefs = length(slimm.references);
    for (uint32_t i=0; i<noOfRefs; ++i)
    {
        if (slimm.references[i].noOfReads == 0)
            continue;
        if (slimm.references[i].covPercent() > slimm.covCutoff &&
            slimm.references[i].noOfReads > slimm.options.minReads)
            slimm.validRefs.insert(i);
    }
    
    uint32_t totalUniqueReads = 0;
    for (auto it= slimm.reads.begin(); it != slimm.reads.end(); ++it)
    {
        it->second.update(slimm.validRefs, slimm.references);
        if(it->second.isUniq(slimm.matchedTaxa, slimm.validRefs))
        {

            __int32 rID = (it->second.targets[0]).rID;
            uint32_t binNo = it->second.targets[0].positions[0];
            
            ++slimm.references[rID].uniqCov2.noOfNonZeroBins;
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
        << current_ref.covDepth() << "\t"
        << current_ref.uniqCovPercent() << "\t"
        << current_ref.uniqCovPercent2() << "\n";
    }
    features_file.close();
}

inline void getReadLCACount(Slimm & slimm,
                            TNodes const & nodes)
{
    for (auto it= slimm.reads.begin(); it != slimm.reads.end(); ++it)
    {
        if (it->second.sumRefLengths > 0)
        {
            std::set<uint32_t> taxaIDs;
            std::set<uint32_t> refIDs;
            
            size_t len = it->second.targets.size();
            for (size_t i=0; i < len; ++i)
            {
                __int32 refID = (it->second.targets[i]).rID;
                taxaIDs.insert(slimm.matchedTaxa[refID]);
                refIDs.insert(refID);
            }
            uint32_t lcaTaxaID = getLCA(taxaIDs, nodes);
            // If taxaID already exists
            if (slimm.taxaID2ReadCount.count(lcaTaxaID) == 1)
                ++slimm.taxaID2ReadCount[lcaTaxaID];
            else   // first time for taxaID
                slimm.taxaID2ReadCount[lcaTaxaID] = 1;
            //add the contributing children references to the taxa
            slimm.taxaID2Children[lcaTaxaID].insert(refIDs.begin(), refIDs.end());
        }
    }
    //add the sum of read counts of children all ancestors of the LCA
    TIntIntMap tID2ReadCountCopy = slimm.taxaID2ReadCount;
    for (auto t2rc : tID2ReadCountCopy)
    {
        uint32_t currentTaxaID = t2rc.first;
        std::set<uint32_t> refIDs = slimm.taxaID2Children[t2rc.first];
        while (nodes.count(currentTaxaID) == 1 && currentTaxaID != 0)
        {
            currentTaxaID = (nodes.at(currentTaxaID)).first;
            if (slimm.taxaID2ReadCount.count(currentTaxaID) == 1)
                ++slimm.taxaID2ReadCount[currentTaxaID];
            else
                slimm.taxaID2ReadCount[currentTaxaID] = 1;
            //add the contributing children references to the taxa
            slimm.taxaID2Children[currentTaxaID].insert(refIDs.begin(), refIDs.end());
        }
    }

}

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
    std::cout<<std::endl<<"Mean = " << m <<" SD = " << sd <<" Cutoff = " << cutoff <<std::endl;
    
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
        abundunceFile   << count << "\t"
                        << candidateName << "\t"
                        << tID.first << "\t"
                        << slimm.taxaID2ReadCount.at(tID.first) << "\t"
                        << tID.second << "\t"
                        << relAbundance << "\n";
        ++count;
    }
    
    // add the remaining  matching reads to unknowns.
    abundunceFile   << count << "\tunknown_"<<rank<< "(multiple)" << "\t0\t"
    << unknownReads << "\t0\t" << unknownAbundance << "\n";
    abundunceFile.close();
    std::cout<< faild_count <<" bellow cutoff ("<< cutoff <<") ...";
}

void getFilesInDirectory(StringList &inputPaths, std::string directory)
{
#ifdef WINDOWS
    HANDLE dir;
    WIN32_FIND_DATA file_data;

    if ((dir = FindFirstFile((directory + "/*").c_str(),
                             &file_data)) == INVALID_HANDLE_VALUE)
        return; /* No files found */

    do
    {
        const std::string file_name = file_data.cFileName;
        const std::string full_file_name = directory + "/" + file_name;
        const bool is_directory =
            (file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;

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
    while ((ent = readdir(dir)) != NULL)
    {
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

#endif /* slimm_h */
