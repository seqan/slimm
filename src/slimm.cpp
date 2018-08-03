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


#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>

#include "timer.hpp"
#include "misc.hpp"
#include "file_helper.hpp"
#include "reference_contig.hpp"
#include "read_stat.hpp"

#include "slimm.hpp"

using namespace seqan;

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------
void setupArgumentParser(ArgumentParser & parser, arg_options const & options)
{
    // Setup ArgumentParser.
    setAppName(parser, "slimm");
    setShortDescription(parser, "Species Level Identification of Microbes from Metagenomes");
    setCategory(parser, "Metagenomics");

    setDateAndVersion(parser);
    setDescription(parser);
    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIDB\\fP\" \"\\fIIN\\fP\"");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "DB"));
    setValidValues(parser, 0, ".sldb");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "IN"));

    // The output file argument.
    addOption(parser, ArgParseOption("o", "output-prefix", "output path prefix.", ArgParseArgument::OUTPUT_PREFIX));

    addOption(parser, ArgParseOption("w", "bin-width", "Set the width of a single bin in neuclotides.",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("mr", "min-reads", "Minimum number of matching reads to consider a reference present.",
                                     ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("r", "rank", "The taxonomic rank of identification", ArgParseOption::STRING));
    setValidValues(parser, "rank", options.rankList);
    setDefaultValue(parser, "rank", options.rank);

    setDefaultValue(parser, "bin-width", options.bin_width);
    setDefaultValue(parser, "min-reads", options.min_reads);

    addOption(parser, ArgParseOption("cc", "cov-cut-off", "the quantile of coverages to use as a cutoff smaller value means bigger threshold.",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));

    setMinValue(parser, "cov-cut-off", "0.0");
    setMaxValue(parser, "cov-cut-off", "1.0");
    setDefaultValue(parser, "cov-cut-off", options.cov_cut_off);

    addOption(parser, ArgParseOption("ac", "abundance-cut-off", "do not report abundances below this value",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    setMinValue(parser, "abundance-cut-off", "0.0");
    setMaxValue(parser, "abundance-cut-off", "10.0");
    setDefaultValue(parser, "abundance-cut-off", options.abundance_cut_off);


    addOption(parser,
              ArgParseOption("d", "directory", "Input is a directory."));
    addOption(parser,
              ArgParseOption("ro", "raw-output", "Output raw reference statstics"));
    addOption(parser,
              ArgParseOption("v", "verbose", "Enable verbose output."));

    // Add Examples Section.
    addTextSection(parser, "Examples");

    addListItem(parser,
                "\\fBslimm\\fP \\fB-o\\fP "
                "\\fIslimm_reports/\\fP \\fIslimm_db_5K.sldb\\fP \\fIexample.bam\\fP",
                "get taxonomic profile from \"\\fIexample.bam\\fP\" "
                "and write it to a tsv file under \"\\fIslimm_reports/\\fP\" directory.");

    addListItem(parser,
                "\\fBslimm\\fP \\fB-d\\fP \\fB-o\\fP "
                "\\fIslimm_reports/\\fP \\fIslimm_db_5K.sldb\\fP \\fIexample-dir/\\fP",
                "get taxonomic profiles from individual SAM/BAM files "
                "located under \"\\fIexample-dir/\\fP\" and write them to tsv files "
                "under \"\\fIslimm_reports/\\fP\" directory with their corsponding file names.");
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------
ArgumentParser::ParseResult
parseCommandLine(ArgumentParser & parser, arg_options & options, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "bin-width"))
        getOptionValue(options.bin_width, parser, "bin-width");

    if (isSet(parser, "min-reads"))
        getOptionValue(options.min_reads, parser, "min-reads");

    if (isSet(parser, "rank"))
        getOptionValue(options.rank, parser, "rank");

    if (isSet(parser, "cov-cut-off"))
        getOptionValue(options.cov_cut_off, parser, "cov-cut-off");

    if (isSet(parser, "abundance-cut-off"))
        getOptionValue(options.abundance_cut_off, parser, "abundance-cut-off");

    if (isSet(parser, "verbose"))
        getOptionValue(options.verbose, parser, "verbose");

    if (isSet(parser, "directory"))
        options.is_directory = true;

    if (isSet(parser, "raw-output"))
        options.raw_output = true;

    getArgumentValue(options.database_path, parser, 0);
    getArgumentValue(options.input_path, parser, 1);

    getOptionValue(options.output_prefix, parser, "output-prefix");
    if (!isSet(parser, "output-prefix"))
        options.output_prefix = options.input_path;

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
    arg_options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(parser, options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    return get_taxonomic_profile(options);
}
