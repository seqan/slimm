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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
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


#include <sstream>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>


using namespace seqan;
typedef seqan::StringSet <Dna5String> SequenceList;
typedef seqan::StringSet <CharString> StringList;

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    bool            verbose;
    CharString      genomesFile;
    CharString      genomesIndexFile;
    __uint32        noOfGenomes;
    __uint64        longestGenomeLen;
    __uint64        concatGenomesLen;
    Options() :
        verbose(false),
        noOfGenomes(),
        longestGenomeLen(),
        concatGenomesLen()
    {}
};

// ----------------------------------------------------------------------------
// Class SlimmIndexer
// ----------------------------------------------------------------------------
template <typename TString, typename TSpec=void>
struct SlimmIndexer
{
    typedef StringSet<CharString>   TGenomeIds;
    typedef StringSet<TString, Owner<ConcatDirect<> > >   TFastaEntry;
    Options const &     options;
    TFastaEntry         genomes;
    TGenomeIds          genomeIds;

    SlimmIndexer(Options const & options) :
        options(options)
    {
    }
};
// --------------------------------------------------------------------------
// Function GetFilesInDirectory() - returns a list of paths given a directory
// --------------------------------------------------------------------------

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


// --------------------------------------------------------------------------
// Function loadReferences()
// --------------------------------------------------------------------------
template <typename TString, typename TSpec=void>
void loadReferences(SlimmIndexer<TString, TSpec> & me)
{
    SeqFileIn genomesFileIn;
    
    if (me.options.verbose)
        std::cerr << "Loading references genomes ... \t\t\t" << std::flush;
    
    if (!open(genomesFileIn, toCString(me.options.genomesFile)))
        throw RuntimeError("Error while opening the reference file(s).");

    try
    {
        readRecords(me.genomeIds, me.genomes, genomesFileIn);
    }
    catch (BadAlloc const & )
    {
        throw RuntimeError("Insufficient memory to load the reference.");
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
    }
    if (me.options.verbose)
        std::cerr << " OK\n";
}



// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser)
{
    // Set short description, version, and date.
    setAppName(parser, "slimm_indexer");
    setShortDescription(parser, "Builds a fasta index for SLIMM application");
    setVersion(parser, "0.1");  
    setDate(parser, "August 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIgenomes.fa genome\\fP>");
    addDescription(parser, "Builds a fasta index for SLIMM application. In "
                   "addition it will collect intrinsic information about each "
                   "genome in the input fasta file.");

    // The input file.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "genome.fa"));
    setValidValues(parser, 0, "fasta fa");
    // The output file argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::OUTPUT_FILE, "genome"));

    addOption(parser, ArgParseOption("v", "verbose", "Displays verbose output."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBslimm_indexer\\fP \\fB-v\\fP \\fIexample.fa\\fP \\fIindices/example\\fP",
    "index the sequences in\"example.fa\" and save the indices with prefix \"indices/example\".");
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    
    if (res != ArgumentParser::PARSE_OK)
        return res;
    // Parse genomes input file.
    getArgumentValue(options.genomesFile, parser, 0);

    // Parse the index prefix.
    getArgumentValue(options.genomesIndexFile, parser, 1);

    // Parse verbose output option.
    getOptionValue(options.verbose, parser, "verbose");

    return ArgumentParser::PARSE_OK;
}


// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.
int main(int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser;
    Options options;
    setupArgumentParser(parser);
    ArgumentParser::ParseResult parsedOptions = parseCommandLine(options, parser, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.

    if (parsedOptions != seqan::ArgumentParser::PARSE_OK)
        return parsedOptions == seqan::ArgumentParser::PARSE_ERROR;

    typedef FMIndexConfig<void, unsigned> TConfig;
    typedef String<Dna5> TString;
    typedef StringSet<TString> TStringSet;
    typedef Index<TStringSet, FMIndex<> > TIndex;
    
    SlimmIndexer<TString> slimmIndexer(options);

    loadReferences(slimmIndexer);
    if (options.verbose)
        std::cerr << "Checking for existing indices ... \t\t" << std::flush;
   
    TIndex fmIndex;
    // Check if there are existing indices built with the same prefix
    if(open(fmIndex, toCString(options.genomesIndexFile)))
    {
        if (options.verbose)
            std::cerr << " OK\n";
        std::cerr << "A previouslly built indices found with the same prefix. "
        "You can:\n 1. Use the existing indices OR\n 2. Use a deferent prefix "
        "OR\n 3. Delete the existing indices.\n" << std::flush;
        return 1;
    }
    
    if (options.verbose)
    {
        std::cerr << "Building the index ... \t\t\t\t\t" << std::flush;
    }
    clear(fmIndex);
    setValue(fmIndex.text, slimmIndexer.genomes);
    std::cout << "No of sequences: " << length(slimmIndexer.genomes) << std::endl;

    try
    {
        // start on demand index construction.
        typename Iterator<TIndex, TopDown<> >::Type it(fmIndex);
        ignoreUnusedVariableWarning(it);
    }
    catch (BadAlloc const &)
    {
        throw RuntimeError("There is no enough memory to index the references.");
        return 1;
    }
    catch (IOError const & )
    {
        throw RuntimeError("There is no enough disk space to build the index");
        return 1;
    }
    if (options.verbose)
        std::cerr << " OK\n";
    if (options.verbose)
        std::cerr << "Saving the index ... \t\t\t\t\t" << std::flush;
    try
    {
//        indexRequire(fmIndex, EsaBwt());       // for (Super-)MaxRepeats iterators
        if (!save(fmIndex, toCString(options.genomesIndexFile)))
            throw RuntimeError("Error while saving the reference index file.");
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    if (options.verbose)
        std::cerr << " OK\n";    
    return 0;
}