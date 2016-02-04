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

#include <seqan/basic.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>


using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

typedef seqan::StringSet <Dna5String> SequenceList;
typedef seqan::StringSet <CharString> StringList;

// Calculates log2 of number.  
float log_2(float n)  
{  
    // log(n)/log(2) is log2.  
    return log(n)/log(2);  
}

// returns a vector after spliting a string in two chunks
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
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

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.
template <typename T>
std::string NumberToString ( T Number )
{
  std::stringstream ss;
  ss<<Number;
  return ss.str();
}

template <typename T>
T StringToNumber ( const std::string &Text )
{
  std::stringstream ss(Text);
  T result;
  return ss >> result ? result : 0;
}


class AlignmentEntry
{
public:
  std::string readId;
  std::string referenceGI;
  unsigned beginPos;
  unsigned int readLength; 
  unsigned int editDistance; 
  seqan::String<CigarElement<> > cigarString; 
  AlignmentEntry() :
  readId(""),
  referenceGI(""),
  beginPos(0),
  readLength(0),
  editDistance(0)
  {}
};

class ReadAlignmnetProperty
{
public:
  unsigned beginPos;
  float score;
  ReadAlignmnetProperty():
  beginPos(0),
  score(0.0)
  {}
};

class TaxaProperty
{
public:
  std::string name;
  unsigned int length; 
  unsigned int neighbors_lv1; 
  unsigned int neighbors_lv2; 
  unsigned int neighbors_lv3; 
  float GCContent;
  TaxaProperty():
  name(""),
  length(0),
  neighbors_lv1(0),
  neighbors_lv2(0),
  neighbors_lv3(0),
  GCContent(50.0)
  {}
};
class ReferenceEntry
{
public:
  std::unordered_map <std::string, std::vector<ReadAlignmnetProperty> > reads;
  unsigned noOfUniqueReads;
  unsigned noOfReads;
  unsigned noOfBasePairs;
  int taxid;
  float validity;
  float mappingError; 
  float homogeneity;
  ReferenceEntry():
  noOfUniqueReads(0),
  // noOfUniqueReadsByTaxid(0),
  noOfReads(0),
  noOfBasePairs(0),
  taxid(0),
  validity(0.0),
  mappingError(0.0), 
  homogeneity(0.0)
  {}
};

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


float CalculateAlignmentScore(seqan::String<CigarElement<> > cigar, int editDistance, unsigned readLen)
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



void GetFilesInDirectory(StringList &input_paths, std::string directory)
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
      if(full_file_name.find(".sam") != std::string::npos || full_file_name.find(".bam") != std::string::npos )
        seqan::appendValue(input_paths, full_file_name);
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


      if(full_file_name.find(".sam") != std::string::npos || full_file_name.find(".bam") != std::string::npos )
        seqan::appendValue(input_paths, full_file_name);
// out.push_back(full_file_name);
    }
    closedir(dir);
#endif
} // GetFilesInDirectory

struct AppOptions
{
// Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
  int verbosity;

  bool isDirectory;
// The first (and only) argument of the program is stored here.
  seqan::CharString input_path;
  seqan::CharString output_path;
  AppOptions() :
  verbosity(1),
  isDirectory(false)
  {}
};


// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
// Setup ArgumentParser.
  seqan::ArgumentParser parser("slimm");
// Set short description, version, and date.
  setShortDescription(parser, "Species Level Identification of Microbes from Metagenomes - Investigates which microbes are present from a SAM alignment file.");
  setVersion(parser, "0.1");
  setDate(parser, "August 2014");

// Define usage line and long description.
  addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIIN\\fP\" \"\\fIOUT\\fP\"");
  addDescription(parser, "Investigates which microbes are present from a BAM/SAM alignment file .");

// The input file/directory argument.
  addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "IN"));
// The output file argument.
  addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));

  addOption(parser, seqan::ArgParseOption("d", "directory", "Input is a directory."));
  addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
  addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));

// Add Examples Section.
  addTextSection(parser, "Examples");
  addListItem(parser, "\\fBslimm\\fP \\fB-v\\fP \\fIexample.bam\\fP \\fIfeatures.csv\\fP",
    "get the microbes present in \"example.bam\" alignment file and write the result to \"features.csv\".");
  addListItem(parser, "\\fBslimm\\fP \\fB-d\\fP \\fIexample-dir\\fP \\fIfeatures.csv\\fP",
    "get the microbes present in the SAM/BAM files located under \"example-dir\" and write the result to \"features.csv\".");

// Parse command line.
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

// Only extract  options if the program will continue after parseCommandLine()
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res;

// Extract option values.
  if (isSet(parser, "quiet"))
    options.verbosity = 0;
  if (isSet(parser, "verbose"))
    options.verbosity = 2;
  if (isSet(parser, "directory"))
    options.isDirectory = true;


  seqan::getArgumentValue(options.input_path, parser, 0);
  seqan::getArgumentValue(options.output_path, parser, 1);

  return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main2(int argc, char const ** argv)
{
// Parse the command line.
  seqan::ArgumentParser parser;
  AppOptions options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

// If there was an error parsing or built-in argument parser functionality
// was triggered then we exit the program.  The return code is 1 if there
// were errors and 0 if there were none.

  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;

// Check if the path exists
  std::ifstream ifile(toCString(options.input_path));
  if (!ifile) {
    std::cout<<options.input_path<<" is not a valid path\n";
    return 1; 
  }

// Prepare the 
  StringList input_paths;
  int numFiles = 1, recCount = 0, badRecCount = 0, fileCount = 0;
  if (options.isDirectory)
  {

    GetFilesInDirectory(input_paths, toCString(options.input_path));
    numFiles = length(input_paths);
    std::cout<<numFiles<<": SAM/BAM Files found under the directory: "<<options.input_path<<"! \n" ;
  }
  else
  {
    if (is_file(toCString(options.input_path)))
      seqan::appendValue(input_paths, toCString(options.input_path));
    else 
    {
      std::cout<<options.input_path<<" is not a file use -d option for a directory.\n";
      return 1;
    }
  }

  Iterator<StringList>::Type it = begin(input_paths);
  Iterator<StringList, Standard>::Type itEnd = end(input_paths); //same Iterator as above

  std::stringstream ss;    
  std::ofstream sam_extract_file;    
  seqan::BamAlignmentRecord record;

  std::unordered_map <unsigned, AlignmentEntry > alignments; //all the alignment information in the sam/bam files
  std::unordered_map <std::string, ReferenceEntry > references; //all the references in the sam/bam files
  std::unordered_map <std::string, std::set<std::string> > readTable; //all the reads in the sam/bam files
  std::vector<std::string> uniqueReads;
  std::vector<std::string> uniqueReadsByTaxid;

  unsigned alignment_index = 0 ;
  for(; it != itEnd; goNext(it))
  {
    badRecCount = 0;
    fileCount ++;
    seqan::CharString curr_file;

    std::cout<<"Reading "<<fileCount<<" of "<< numFiles<<" files ...";
    // read the original file
    seqan::BamFileIn bamIn(toCString(value(it)));
    seqan::BamHeader bamHead;
    readHeader(bamHead, bamIn);

    // Read alignment records from the sam/bam files
    while(!atEnd(bamIn))
    {
      try {
        readRecord(record, bamIn);
        if(!(hasFlagUnmapped(record)))
        {
          seqan::CharString contigName = getContigName(record, bamIn);
          seqan::StringSet<CharString> splitedContigName;
          strSplit(splitedContigName, contigName, EqualsChar<'|'>(), true, 4);
          // std::cout<< splitedContigName[1] <<std::endl; 
          AlignmentEntry ali_entry;
          ali_entry.readId = toCString(record.qName);
          ali_entry.referenceGI = toCString(splitedContigName[1]);
          ali_entry.beginPos = record.beginPos;
          ali_entry.readLength = length(record.seq);
          ali_entry.cigarString = record.cigar;
          BamTagsDict recordTags(record.tags);
          unsigned id = 0;
          findTagKey(id, recordTags, "NM");
          extractTagValue(ali_entry.editDistance, recordTags, id);

          alignments[alignment_index] = ali_entry;
          alignment_index ++;
          recCount ++;
        }
      }
      catch (Exception const & e)
      {
          std::cout << "ERROR: " << e.what() << std::endl;
          badRecCount ++;
          std::cerr<<"Can't read record!"<<std::endl;
          // If invalid records keep comming and we reach the treshold 
          // stop reading the file, and move on to the next step
          if (badRecCount > 20)
          {
            std::cerr<<"Too many bad records! Moving on to the next file"<<std::endl;
            break;
          }
      }
    }
    std::cout<<"\tDone!" <<std::endl;
  }

  // ============get the gi2taxid mapping ================
  std::cout<<"Loding gi2taxid mapping ...";
  std::unordered_map <unsigned, unsigned> gi2taxid;
  std::ifstream giMap("gi2taxid.map");
  int gi, taxid;
  while (giMap >> gi >> taxid)
  {
    gi2taxid[gi]=taxid;
  }
  giMap.close();
  std::cout<<"\tDone!" <<std::endl;  
  

  std::cout<<"Analysing alignments, reads and references ...";
  for (unsigned int i=0; i< alignments.size(); i++)
  {
    ReadAlignmnetProperty readAliPropertiy;
    readAliPropertiy.beginPos = alignments[i].beginPos; 
    // extract the bigin position of the alignment
    readAliPropertiy.score = CalculateAlignmentScore(alignments[i].cigarString, 
                          alignments[i].editDistance, alignments[i].readLength);

    std::string current_refGI = alignments[i].referenceGI;
    std::string current_readId = alignments[i].readId;

    readTable[current_readId].insert(current_refGI);

    if (references.count(current_refGI) == 0)
    {
      ReferenceEntry refEntry;
      refEntry.reads[current_readId].push_back(readAliPropertiy) ;
      refEntry.mappingError = readAliPropertiy.score;
      refEntry.noOfBasePairs = alignments[i].readLength;
      references[current_refGI] = refEntry;
    }
    else
    {
      references[current_refGI].reads[current_readId].push_back(readAliPropertiy);
      references[current_refGI].mappingError += readAliPropertiy.score;
      references[current_refGI].noOfBasePairs += alignments[i].readLength;
    }
  }
  std::cout<<"\tDone!" <<std::endl;


  unsigned current = 0;
  unsigned total = references.size();

  std::cout<<"computing features of each reference genome ... " <<std::endl;
  typedef std::unordered_map <std::string, ReferenceEntry >::iterator ReferenceIterator;
  for (ReferenceIterator i = references.begin(); i != references.end(); i++)
  {
    current ++;
    std::cout<<current <<" of "<<total<<"\r";
    i->second.taxid = gi2taxid[atoi(toCString(i->first))];
    i->second.noOfReads = i->second.reads.size();
    i->second.mappingError = i->second.mappingError/i->second.reads.size();
  }
  std::cout<<"\t\t\t\tDone!" <<std::endl;

  std::cout<<"Calculating Unique Reads ...";
  unsigned totalUniqueReads = 0;
  typedef std::unordered_map <std::string, std::set<std::string> >::iterator ReadIterator;
  for (ReadIterator i = readTable.begin(); i != readTable.end(); i++)
  {
    std::set<std::string>::iterator it = i->second.begin();
    if (i->second.size() == 1)
    {
      references[*it].noOfUniqueReads++;
      totalUniqueReads ++;
    }

  }
  std::cout<<"\tDone!" <<std::endl;
  std::cout<<"Total number of Unique reads = "<< totalUniqueReads <<std::endl;

  // ============get the taxon properties ================
  std::cout<<"Loding additional taxon related features ...";
  std::unordered_map <int, TaxaProperty> taxonProperties;
  std::ifstream taxidMap("taxidFeatures.map");
  std::string line;
  while (std::getline(taxidMap, line))
  {
    std::vector<std::string> chunks = split(line, '\t');
    TaxaProperty tp;
    tp.name = chunks[1];
    tp.length = atof(toCString(chunks[2]));
    tp.GCContent = atof(toCString(chunks[3]));
    tp.neighbors_lv1 = atoi(toCString(chunks[4]));
    tp.neighbors_lv2 = atoi(toCString(chunks[5]));
    tp.neighbors_lv3 = atoi(toCString(chunks[6]));

    taxonProperties[atoi(toCString(chunks[0]))]=tp;
  }
  taxidMap.close();
  std::cout<<"\tDone!" <<std::endl;

  std::cout<<"Writing features to a file ..." <<std::endl;
  std::ofstream features_file;
  features_file.open(toCString(options.output_path));
  
  features_file<<"No.\tCandidate Name\tTaxid\tGen. Length\t"
  "Num. Reads\tUnique Reads\tCoverage\t"
  "Mapping Error\tRel. abundance.\tUnique Ratio\t"
  "GC content\tneighbors_lv1\tneighbors_lv2\tneighbors_lv3\n";

  current = 0;
  for (ReferenceIterator i = references.begin(); i != references.end(); i++)
  {
    current ++;

    int currentTaxid = i->second.taxid;
    std::cout<<current <<" of "<<total<<"\r";

    std::string candidateName = i->first;
    if (taxonProperties[currentTaxid].name.length() > 0)
      candidateName = taxonProperties[currentTaxid].name;

    features_file<<current<<"\t"<<candidateName<<"\t"<<i->second.taxid<<"\t" 
    <<taxonProperties[currentTaxid].length<<"\t"<<i->second.noOfReads<<"\t" 
    <<i->second.noOfUniqueReads<<"\t"
    <<float(i->second.noOfBasePairs)/taxonProperties[currentTaxid].length<<"\t" 
    <<i->second.mappingError<<"\t"
    <<float(i->second.noOfReads)/readTable.size()<<"\t"
    <<float(i->second.noOfUniqueReads)/totalUniqueReads<<"\t"
    <<taxonProperties[currentTaxid].GCContent<<"\t"
    <<taxonProperties[currentTaxid].neighbors_lv1<<"\t"
    <<taxonProperties[currentTaxid].neighbors_lv2<<"\t"
    <<taxonProperties[currentTaxid].neighbors_lv3<<"\n";
  }
  features_file.close();
  std::cout<<"\t\t\t\tDone!" <<std::endl;


  AlignmentEntry sampleAlignment = alignments[30454];
  std::cout<<"Map size: "<<alignments.size()<<std::endl;
  std::cout<<"Number of references: "<<references.size()<<std::endl;

  std::cout<<recCount<<" SAM/BAM alignment records are proccessed."<<std::endl;
  std::cout<<"extracted features are written to: " <<options.output_path<<std::endl;
return 0;
}