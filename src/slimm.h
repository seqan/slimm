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


#include "timer.h"
#include "load_mapping_files.h"


#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <sys/stat.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <numeric>
#include <unordered_map>




using namespace seqan;

typedef StringSet <Dna5String> SequenceList;
typedef StringSet <CharString> StringList;
typedef std::unordered_map <uint32_t, std::pair<uint32_t, std::string> > TNodes;


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
    typedef std::vector<std::string>            TList;
    
    TList               rankList        = {
        "all",
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "superkingdom"
    };

    float               covCutOff;
    __intSizeBinWidth   binWidth;
    __intSizeBinWidth   minReads;
    bool                verbose;
    bool                isDirectory;
    bool                outputRaw;
    std::string         rank;
    std::string         inputPath;
    std::string         outputPrefix;
    std::string         mappingDir;
    
    AppOptions() :
    covCutOff(0.99),
    binWidth(0),
    minReads(100),
    verbose(false),
    isDirectory(false),
    outputRaw(false),
    rank("all"),
    inputPath(""),
    outputPrefix(""),
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
    int32_t _noOfNonZeroBins = -1;
public:
    uint32_t                    binWidth;
    uint32_t                    noOfBins;
    std::vector <uint32_t>      binsHeight;
    
    Coverage():
    binWidth(0),
    noOfBins(0)
    {}
    
    Coverage(uint32_t totalLen, uint32_t width)
    {
        binWidth = width;
        noOfBins = totalLen/width + 1;
        binsHeight.resize(noOfBins);
    }
    uint32_t noOfNonZeroBins()
    {
        if (_noOfNonZeroBins == -1)
        {
            _noOfNonZeroBins = (noOfBins - std::count (binsHeight.begin(), binsHeight.end(), 0));
        }
        return _noOfNonZeroBins;
    }

};

class ReferenceContig
{
    float _covPercent = -1;
    float _uniqCovPercent = -1;
    float _uniqCovPercent2 = -1;
    float _covDepth = -1;
    float _uniqCovDepth = -1;
    float _uniqCovDepth2 = -1;
public:
    CharString          refName;
    bool                isValid;
    uint32_t            length;
    uint32_t            noOfReads;
    uint32_t            noOfUniqReads;
    uint32_t            noOfUniqReads2;
    Coverage            cov;
    Coverage            uniqCov;
    Coverage            uniqCov2;
    uint32_t            taxaID;
    float               relAbundance;
    float               relAbundanceUniq;
    float               relAbundanceUniq2;
    
    ReferenceContig():
    refName(""),
    isValid(false),
    length(0),
    noOfReads(0),
    noOfUniqReads(0),
    noOfUniqReads2(0),
    cov(),
    uniqCov(),
    uniqCov2(),
    taxaID(0),
    relAbundance(0.0),
    relAbundanceUniq(0.0),
    relAbundanceUniq2(0.0)
    {}

    //Member functions
    float               covPercent();
    float               uniqCovPercent();
    float               uniqCovPercent2();
    float               covDepth();
    float               uniqCovDepth();
    float               uniqCovDepth2();
};


class TargetRef
{
public:
    __int32                    rID;
    std::vector<uint32_t>      positions;

    //constructer takes a ref id and a position for the first time
    TargetRef(__int32 ref, uint32_t pos)
    {
        rID = ref;
        positions.reserve(10);
        positions.push_back(pos);
    }
};


class Read
{
public:
    std::vector<TargetRef>          targets;
    uint32_t                        sumRefLengths = 0;
    uint32_t                        len = 0;
    
    //checks if all the match points are in the same sequence
    bool isUniq();
    
    
    //checks if all the match points are in the same taxaID
    bool isUniq(std::vector<uint32_t> const & taxaIDs);
    
    void update(std::vector<uint32_t> const & taxaIDs,
                std::set<uint32_t> const & valtaxaIDs,
                std::vector<ReferenceContig> const & references );
    
    void addTarget(int32_t rID, uint32_t binNo);

};

// ----------------------------------------------------------------------------
// Class Slimm
// ----------------------------------------------------------------------------
class Slimm
{
    float _covCutoff = 0.0;
    float _uniqCovCutoff = 0.0;
    int32_t _minUniqReads = -1;
    int32_t _minReads = -1;
public:
    std::vector<ReferenceContig>    references;
    std::vector<uint32_t>           matchedTaxa;
    std::set<uint32_t>              validRefTaxonIDs;
    std::unordered_map<std::string, Read>      reads;
    __intSizeQLength                avgQLength = 0;
    __intSizeGLength                matchedRefsLen = 0;
    __intSizeMatchCount             noOfRefs = 0;
    __intSizeMatchCount             failedByMinRead = 0;
    __intSizeMatchCount             failedByMinUniqRead = 0;
    __intSizeMatchCount             failedByCov = 0;
    __intSizeMatchCount             failedByUniqCov = 0;
    __intSizeQCount                 hitCount = 0;
    __intSizeQCount                 uniqHitCount = 0;
    __intSizeQCount                 noOfMatched = 0;
    __intSizeQCount                 noOfUniqlyMatched = 0;
    __intSizeQCount                 noOfUniqlyMatched2 = 0;
    
    AppOptions                      options;
    
    TIntIntMap                      taxaID2ReadCount;
    TIntFloatMap                    taxaID2Abundance;
    TIntSetMap                      taxaID2Children;
    float                           expCov() const;

    float                           covCutoff();
    float                           uniqCovCutoff();
    __intSizeQCount                 minReads();
    __intSizeQCount                 minUniqReads();
    
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
bool Read::isUniq(std::vector<uint32_t> const & taxaIDs)
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
            refTaxaIDs.insert(refID);
        }
        if (refTaxaIDs.size() > 1)
        {
            return false;
        }
    }
    return true;
}


void Read::update(std::vector<uint32_t> const & taxaIDs,
                  std::set<uint32_t> const & valtaxaIDs,
                  std::vector<ReferenceContig> const & references )
{
    size_t len = targets.size();
    if (len == 0)
        return;
    else
    {
        std::vector<TargetRef> newTargets;
        for (size_t i=0; i < len; ++i)
        {
            uint32_t refID = taxaIDs[(targets[i]).rID];
            if(valtaxaIDs.find(refID) != valtaxaIDs.end())
                newTargets.push_back(targets[i]);
            else
                sumRefLengths -= references[(targets[i]).rID].length;
        }
        targets = newTargets;
    }
}

void Read::addTarget(int32_t rID, uint32_t binNo)
{
    if (targets.empty())
    {
        targets.push_back(TargetRef(rID, binNo));
        return;
    }
    else
    {
        for (auto tar : targets)
        {
            if(tar.rID == rID)
            {
                tar.positions.push_back(binNo);
                return;
            }
        }
        targets.push_back(TargetRef(rID, binNo));
    }
}


// ==========================================================================
// Functions
// ==========================================================================

template <typename Type>
Type getCutoffByQuantile (std::vector<Type> v,float q)
{
    Type total = std::accumulate(v.begin(), v.end(), (Type)0);
    Type cutoff = (Type)0, subTotal = (Type)0;
    std::vector<Type> vals = v;
    std::sort (vals.begin(), vals.end());
    
    uint32_t i= vals.size() - 1;
    while((float(subTotal)/total) < q && i > (Type)0)
    {
        subTotal += vals[i];
        --i;
    }
    cutoff = vals[i];
    
    return cutoff;
}


template <typename Type>
bool greaterThan (Type i,Type j)
{
    return (i>j);
}
template <typename Type>
bool lessThan (Type i,Type j)
{
    return (i<j);
}


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
    T vSum = std::accumulate(v.begin(), v.end(), (T)0);
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

std::string getTSVFileName (const std::string& oPrefix, const std::string& inpfName)
{
    std::string result = getDirectory(oPrefix);
    result.append("/");
    if (length(getFilename(oPrefix)) > 0)
        result.append(getFilename(oPrefix));
    else
        result.append(getFilename(inpfName));
    if ((result.find(".sam") != std::string::npos && result.find(".sam") == result.find_last_of(".")) ||
       (result.find(".bam") != std::string::npos && result.find(".bam")  == result.find_last_of(".")))
        result.replace((result.find_last_of(".")), 4, "");
    result.append(".tsv");
    return result;
}

std::string getTSVFileName (const std::string& oPrefix, const std::string& inpfName, const std::string& rank)
{
    std::string fName = getTSVFileName(oPrefix, inpfName);
    std::string sfx = "_" + rank + "_reported";
    return fName.insert(fName.size()-4, sfx);
}

// ----------------------------------------------------------------------------
// Function setDateAndVersion()
// ----------------------------------------------------------------------------

void setDateAndVersion(ArgumentParser & parser)
{
    setDate(parser, __DATE__);
    setShortDescription(parser, "Species Level Identification of Microbes from Metagenomes");
    setCategory(parser, "Metagenomics");
#if defined(SEQAN_APP_VERSION) && defined(SEQAN_REVISION)
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
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
    return float(cov.noOfNonZeroBins())/cov.noOfBins;
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
inline float getCovDepth(Coverage c)
{
    if(c.noOfNonZeroBins() == 0)
        return 0.0;
    
    std::vector <uint32_t>::iterator it;
    
    //copy the coverage height.
    std::vector <float>  bHeights;

    bHeights.reserve(c.noOfBins);
    for (unsigned int i=0; i<c.noOfBins; ++i )
    {
        bHeights.push_back(float(c.binsHeight[i]));
    }
    return mean(bHeights);
}

float ReferenceContig::covDepth()
{
    if (_covDepth == -1)
        _covDepth = getCovDepth(cov);
    return _covDepth;
}
float ReferenceContig::uniqCovDepth()
{
    if (_uniqCovDepth == -1)
        _uniqCovDepth = getCovDepth(uniqCov);
    return _uniqCovDepth;
}
float ReferenceContig::uniqCovDepth2()
{
    if (_uniqCovDepth2 == -1)
        _uniqCovDepth2 = getCovDepth(uniqCov2);
    return _uniqCovDepth2;
}
float Slimm::covCutoff()
{
    if (_covCutoff == 0.0 && options.covCutOff < 1.0)
    {
        std::vector<float> covs = {};
        covs.reserve(length(references));
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].noOfUniqReads > 0)
            {
                covs.push_back(references[i].covPercent());
            }
        }
        _covCutoff = getCutoffByQuantile<float>(covs, options.covCutOff);
    }
    return _covCutoff;
}
float Slimm::uniqCovCutoff()
{
    if (_uniqCovCutoff == 0.0 && options.covCutOff < 1.0)
    {
        std::vector<float> covs = {};
        covs.reserve(length(references));
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].noOfUniqReads > 0)
            {
                covs.push_back(references[i].uniqCovPercent());
            }
        }
        _uniqCovCutoff = getCutoffByQuantile<float>(covs, options.covCutOff);
    }
    return _uniqCovCutoff;
}
float Slimm::expCov() const
{
    return float(avgQLength * noOfMatched) / matchedRefsLen;
}
__intSizeQCount Slimm::minReads()
{
    if (options.covCutOff == 1.0)
        _minReads = 0;
    if (_minReads == -1)
    {
        std::vector<int> counts = {};
        counts.reserve(length(references));
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].noOfReads > 0)
            {
                counts.push_back(references[i].noOfReads);
            }
        }
        _minReads = getCutoffByQuantile(counts, options.covCutOff);
    }
    return _minReads;
}
__intSizeQCount Slimm::minUniqReads()
{
    if (options.covCutOff == 1.0)
        _minUniqReads = 0;
    if (_minUniqReads == -1)
    {
        std::vector<int> uniqCounts = {};
        uniqCounts.reserve(length(references));
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].noOfUniqReads > 0)
            {
                uniqCounts.push_back(references[i].noOfUniqReads);
            }
        }
        _minUniqReads = getCutoffByQuantile(uniqCounts, options.covCutOff);
    }
    return _minUniqReads;
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
        std::set<uint32_t>::iterator it = parents.begin();
        uint32_t t1 = *it;
        ++it;
        uint32_t t2 = *it;
        bool found = false;
        while (nodes.count(t1) == 1 && t1 != 0)
        {
            t2 = *it;
            while (nodes.count(t2) == 1 && t2 != 0)
            {
                if (t1 == t2)
                {
                    found = true;
                    break;
                }
                t2 = nodes.at(t2).first;
            }
            if (found)
            {
                break;
            }
            t1 = nodes.at(t1).first;
        }
        if (found)
        {
            parents.erase(parents.begin(), std::next(it));
            parents.insert(t1);
        }
        else
        {
            return 0;
        }
    }
    return *(parents.begin());
}

bool getTaxaId(uint32_t &idPosition, CharString refName, std::string idType)
{
    StringList chunks;
    strSplit(chunks, refName, EqualsChar<'|'>());
    //check for slimm taxid
    for (uint32_t i = 0; i <  length(chunks) ; ++i)
    {
        if (chunks[i] == idType)
        {
            idPosition = i + 1;
            return true;
        }
    }
    return false;
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
    if (s.size() == 1)
        return taxaIDs[0];
    else
        return getLCA(s, s, nodes);
}


inline void analyzeAlignments(Slimm & slimm,
                              BamFileIn & bamFile)
{

    uint32_t queryLen = 0;
    BamAlignmentRecord record;
    while (!atEnd(bamFile))
    {
        readRecord(record, bamFile);
        if (hasFlagUnmapped(record) || record.rID == BamAlignmentRecord::INVALID_REFID)
            continue;  // Skip these records.

        uint32_t newQueryLen = length(record.seq);
        queryLen = (newQueryLen == 0) ? queryLen : newQueryLen;
        uint32_t center_position = record.beginPos + (queryLen/2);
        if (center_position > slimm.references[record.rID].length)
        {
            center_position = slimm.references[record.rID].length;
        }
        uint32_t relativeBinNo = center_position/slimm.options.binWidth;

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
    
    if (slimm.hitCount != 0)
    {
        __intSizeGLength concatQLength = 0;
        for (auto it= slimm.reads.begin(); it != slimm.reads.end(); ++it)
        {
            concatQLength += it->second.len;
            if(it->second.isUniq(slimm.matchedTaxa))
            {
                __int32 rID = it->second.targets[0].rID;
                it->second.sumRefLengths += slimm.references[rID].length;
                ++slimm.noOfUniqlyMatched;

                size_t pos_count = (it->second.targets[0]).positions.size();
                slimm.references[rID].noOfReads += pos_count;
                it->second.sumRefLengths += slimm.references[rID].length;
                for (size_t j=0; j < pos_count; ++j)
                {
                    uint32_t binNo = (it->second.targets[0]).positions[j];
                    ++slimm.references[rID].cov.binsHeight[binNo];
                }
                slimm.references[rID].noOfUniqReads += 1;
                slimm.uniqHitCount += 1;
                ++slimm.references[rID].uniqCov.binsHeight[(it->second.targets[0]).positions[0]];
            }
            else
            {
                size_t len = it->second.targets.size();
                for (size_t i=0; i < len; ++i)
                {

                    __int32 rID = it->second.targets[i].rID;
                    it->second.sumRefLengths += slimm.references[rID].length;

                    // ***** all of the matches in multiple pos will be counted *****
                    slimm.references[rID].noOfReads += (it->second.targets[i]).positions.size();
                    for (auto binNo : (it->second.targets[i]).positions)
                    {
                        ++slimm.references[rID].cov.binsHeight[binNo];
                    }
                }
            }
        }
        slimm.noOfMatched = slimm.reads.size();

        slimm.avgQLength = concatQLength/slimm.noOfMatched;
        float totalAb = 0.0;
        for (uint32_t i=0; i<length(slimm.references); ++i)
        {
            if (slimm.references[i].noOfReads > 0)
            {
                ++slimm.noOfRefs;
                slimm.matchedRefsLen += slimm.references[i].length;
                slimm.references[i].relAbundance = float(slimm.references[i].noOfReads * 100)/slimm.hitCount;
                totalAb += slimm.references[i].relAbundance/slimm.references[i].length;
            }
            else
            {
                slimm.references[i].relAbundance = 0.0;
            }
        }
        for (uint32_t i=0; i<length(slimm.references); ++i)
        {
            if (slimm.references[i].noOfReads > 0)
            {
                slimm.references[i].relAbundance = (slimm.references[i].relAbundance * 100) / (totalAb*slimm.references[i].length);
            }
        }
        
        totalAb = 0.0;
        for (uint32_t i=0; i<length(slimm.references); ++i)
        {
            if (slimm.references[i].noOfUniqReads > 0)
            {
                slimm.references[i].relAbundanceUniq = float(slimm.references[i].noOfUniqReads * 100)/slimm.uniqHitCount;
                totalAb += slimm.references[i].relAbundanceUniq/slimm.references[i].length;
            }
            else
            {
                slimm.references[i].relAbundanceUniq = 0.0;
            }
        }
        for (uint32_t i=0; i<length(slimm.references); ++i)
        {
            if (slimm.references[i].noOfUniqReads > 0)
            {
                slimm.references[i].relAbundanceUniq = (slimm.references[i].relAbundanceUniq * 100) / (totalAb*slimm.references[i].length);
            }            
        }
    }
}

inline void filterAlignments(Slimm & slimm)
{
    uint32_t noOfRefs = length(slimm.references);
    for (uint32_t i=0; i<noOfRefs; ++i)
    {
        
        if (slimm.references[i].noOfReads == 0)
            continue;
        if (
            slimm.references[i].covPercent() >= slimm.covCutoff() &&
            slimm.references[i].uniqCovPercent() >= slimm.uniqCovCutoff() &&
            true
            )
            slimm.validRefTaxonIDs.insert(slimm.matchedTaxa[i]);
        else
        {
            if(slimm.references[i].uniqCovPercent() < slimm.uniqCovCutoff())
            {
                ++slimm.failedByUniqCov;
            }
            if(slimm.references[i].noOfReads < slimm.options.minReads)
            {
                ++slimm.failedByMinRead;
            }
            if(slimm.references[i].covPercent() < slimm.covCutoff())
            {
                ++slimm.failedByCov;
            }
        }

    }
    
    for (auto it= slimm.reads.begin(); it != slimm.reads.end(); ++it)
    {
        it->second.update(slimm.matchedTaxa, slimm.validRefTaxonIDs, slimm.references);
        if(it->second.isUniq(slimm.matchedTaxa))
        {
            __int32 rID = (it->second.targets[0]).rID;
            slimm.references[rID].noOfUniqReads2 += 1;
            slimm.noOfUniqlyMatched2 += 1;
            uint32_t binNo = (it->second.targets[0]).positions[0];
            ++slimm.references[rID].uniqCov2.binsHeight[binNo];
        }
    }
}


inline void writeToFile(std::string & filePath,
                        std::vector<ReferenceContig> & refList,
                        TIntStrMap const & taxaID2name)
{
    std::ofstream features_file;
    features_file.open(filePath);
    
    features_file <<    "No.\t"
    "CandidateName\t"
    "Taxid\t"
    "NoOfReads\t"
    "RelAbundance\t"
    "RelAbundanceUniq\t"
    "RelAbundanceUniq2\t"
    "GenomeLength\t"
    "NoOfUniqueReads\t"
    "NoOfUniqueReads2\t"
    
    "NoOfBins\t"
    "noOfNonZeroBins\t"
    "noOfNonZeroBinsUniq\t"
    "noOfNonZeroBinsUniq2\t"
    
    
    "CoverageDepth\t"
    
    "UniqCoverageDepth\t"
    "UniqCoverageDepth2\t"
    
    "MappingError\t"
    "CoveragePercentage\t"
    "UniqueCoveragePercentage\t"
    "UniqueCoveragePercentage2\n";
    
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
        << current_ref.noOfReads << "\t"
        << current_ref.relAbundance << "\t"
        << current_ref.relAbundanceUniq << "\t"
        << current_ref.relAbundanceUniq2 << "\t"
        << current_ref.length << "\t"
        << current_ref.noOfUniqReads << "\t"
        << current_ref.noOfUniqReads2 << "\t"
        
        << current_ref.cov.noOfBins << "\t"
        << current_ref.cov.noOfNonZeroBins() << "\t"
        << current_ref.uniqCov.noOfNonZeroBins() << "\t"
        << current_ref.uniqCov2.noOfNonZeroBins() << "\t"
        
        
        << current_ref.covDepth() << "\t"
        
        << current_ref.uniqCovDepth() << "\t"
        << current_ref.uniqCovDepth2() << "\t"
        
        << "NA"<< "\t"
        << current_ref.covPercent() << "\t"
        << current_ref.uniqCovPercent() << "\t"
        << current_ref.uniqCovPercent2() << "\n";
    }
    features_file.close();
}


inline void getReadLCACount(Slimm & slimm,
                            TNodes const & nodes)
{
    // put the non-unique read to upper taxa.
    for (auto it= slimm.reads.begin(); it != slimm.reads.end(); ++it)
    {
        if(!(it->second.isUniq(slimm.matchedTaxa)))
        {
            uint32_t lcaTaxaID = 0;
            std::set<uint32_t> refIDs = {};
            std::set<uint32_t> taxaIDs;
            size_t len = it->second.targets.size();
            for (size_t i=0; i < len; ++i)
            {
                __int32 refID = (it->second.targets[i]).rID;
                taxaIDs.insert(slimm.matchedTaxa[refID]);
                refIDs.insert(refID);
            }
            lcaTaxaID = getLCA(taxaIDs, nodes);
            slimm.taxaID2Children[lcaTaxaID].insert(refIDs.begin(), refIDs.end());
            // If taxaID already exists
            if (slimm.taxaID2ReadCount.count(lcaTaxaID) == 1)
                ++slimm.taxaID2ReadCount[lcaTaxaID];
            else   // first time for taxaID
                slimm.taxaID2ReadCount[lcaTaxaID] = 1;
            //add the contributing children references to the taxa
        }
    }
    //add the sum of read counts of children all ancestors of the LCA
    TIntIntMap tID2ReadCountCopy = slimm.taxaID2ReadCount;
    for (auto t2rc : tID2ReadCountCopy)
    {
        uint32_t currentTaxaID = t2rc.first;
        uint32_t readCount = tID2ReadCountCopy[currentTaxaID];
        std::set<uint32_t> refIDs = slimm.taxaID2Children[currentTaxaID];
        while (nodes.count(currentTaxaID) == 1 && currentTaxaID != 0)
        {
            currentTaxaID = (nodes.at(currentTaxaID)).first;
            if (slimm.taxaID2ReadCount.count(currentTaxaID) >= 1)
                slimm.taxaID2ReadCount[currentTaxaID] += readCount;
            else
                slimm.taxaID2ReadCount[currentTaxaID] = readCount;
            //add the contributing children references to the taxa
            slimm.taxaID2Children[currentTaxaID].insert(refIDs.begin(), refIDs.end());
        }
    }
    
    
    for (uint32_t i=0; i<length(slimm.references); ++i)
    {

        if (slimm.references[i].noOfUniqReads2 > 0)
        {
            uint32_t currentTaxaID = slimm.references[i].taxaID;
            uint32_t uniqCount = slimm.references[i].noOfUniqReads2;
            while (nodes.count(currentTaxaID) == 1 && currentTaxaID != 0)
            {
                if (slimm.taxaID2ReadCount.count(currentTaxaID) >= 1)
                    slimm.taxaID2ReadCount[currentTaxaID] += uniqCount;
                else
                    slimm.taxaID2ReadCount[currentTaxaID] = uniqCount;
                slimm.taxaID2Children[currentTaxaID].insert(i);
                //add the contributing children references to the taxa
                currentTaxaID = (nodes.at(currentTaxaID)).first;
            }
        }
    }
}


inline void writeAbundance(Slimm & slimm,
                           TNodes & nodes,
                           TIntStrMap const & taxaID2name,
                           uint32_t const & noReadsAtRank,
                           std::string const & currFile,
                           std::string const & rank)
{

    std::string tsvFile = getTSVFileName(toCString(slimm.options.outputPrefix), currFile, rank);
    std::ofstream abundunceFile;
    abundunceFile.open(tsvFile);

    uint32_t unknownReads = slimm.noOfMatched-noReadsAtRank;
    uint32_t count = 0;
    uint32_t faild_count = 0;
    float faildAbundunce = 0.0;
    TIntFloatMap cladeCov;
    TIntFloatMap cladeAbundance;
    count = 1;
    abundunceFile<<"No.\tName\tTaxid\tNoOfReads\tRelativeAbundance\tContributers\tCoverage\n";
    for (auto tID : slimm.taxaID2ReadCount) {
        if (rank == nodes[tID.first].second)
        {
            uint32_t cLength = 0;
            uint32_t noOfContribs = 0;
            std::set<uint32_t>::iterator it;
            for (it=slimm.taxaID2Children.at(tID.first).begin();
                 it!=slimm.taxaID2Children.at(tID.first).end(); ++it)
            {
                cLength += slimm.references[*it].length;
                ++noOfContribs;
            }
            cLength = cLength/noOfContribs;
            float cov = float(tID.second * slimm.avgQLength)/cLength;

            float relAbundance = float(tID.second)/(slimm.noOfMatched) * 100;
            std::unordered_map <uint32_t, std::string>::const_iterator it2 =
            taxaID2name.find (tID.first);
            if (relAbundance == 0.001 || cov < slimm.covCutoff())
            {
                faildAbundunce += relAbundance;
                ++faild_count;
                continue;
            }
            seqan::CharString candidateName = "Organism name not found";
            if (it2 != taxaID2name.end())
                candidateName = (taxaID2name.at(tID.first));
            abundunceFile   << count << "\t"
            << candidateName << "\t"
            << tID.first << "\t"
            << tID.second << "\t"
            << relAbundance << "\t"
            << slimm.taxaID2Children.at(tID.first).size() << "\t"
            << cov << "\n";
            ++count;
        }

    }
    float unknownAbundance = float(unknownReads)/ slimm.noOfMatched * 100 + faildAbundunce;

    abundunceFile   << count << "\tunknown_"<< rank<< "(multiple)" << "\t0\t"
    << unknownReads << "\t" << unknownAbundance << "\t0\t0.0\n";
    abundunceFile.close();
    std::cout<< std::setw (15) << rank <<" level: "<< faild_count <<" bellow cutoff ("<< 0.001 <<") ...\n";
}


inline void writeAbundance(Slimm & slimm, TNodes & nodes, TIntStrMap const & taxaID2name, std::string const & currFile)
{
    std::vector<std::string> const allRanks= {"species", "genus", "family", "order", "class", "phylum", "superkingdom"};
    std::vector<std::string> consideredRanks = {slimm.options.rank};
    if(slimm.options.rank == "all")
    {
        consideredRanks.resize(7);
        consideredRanks = allRanks;
    }

    // calculate the total number of reads matching uniquily at that species level.
    std::map<std::string, uint32_t>  noReadsAtRank = {{"all", 0}, {"species", 0}, {"genus", 0}, {"family", 0}, {"order", 0}, {"class", 0}, {"phylum", 0}, {"superkingdom", 0}};
    for (auto tID : slimm.taxaID2ReadCount)
    {
        if (noReadsAtRank.find(nodes[tID.first].second) != noReadsAtRank.end() )
        {
            noReadsAtRank[nodes[tID.first].second] +=  tID.second ;     // found
        }
    }

    for (std::string rank : consideredRanks)
    {
        writeAbundance(slimm, nodes, taxaID2name, noReadsAtRank[rank], currFile, rank);
    }
}
void getFilesInDirectory(std::vector<std::string> & inputPaths, std::string directory)
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
           inputPaths.push_back(full_file_name);
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
            inputPaths.push_back(full_file_name);
        //appendValue(inputPaths, full_file_name);
    }
    closedir(dir);
#endif
} // getFilesInDirectory

#endif /* slimm_h */
