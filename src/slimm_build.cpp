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
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>


#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include "misc.hpp"

using namespace seqan;

// ----------------------------------------------------------------------------
// Class arg_options
// ----------------------------------------------------------------------------
struct arg_options
{
    uint32_t                     batch;
    bool                         verbose;
    std::string                  fasta_path;
    std::string                  nodes_path;
    std::string                  names_path;
    std::string                  output_path;
    std::vector<std::string>     ac__taxid_paths;

    arg_options() : batch(1000000),
                    verbose(false),
                    fasta_path(),
                    nodes_path(),
                    names_path(),
                    output_path("slimm_db.sldb"),
                    ac__taxid_paths() {}
};

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------
void setupArgumentParser(ArgumentParser & parser, arg_options const & options)
{
    // Setup ArgumentParser.
    setAppName(parser, "slimm-build");
    setShortDescription(parser, "gets a reduced taxonomic information given a multi-fasta file using accession numbers");
    setCategory(parser, "Metagenomics");

    setDateAndVersion(parser);
    setDescription(parser);
    // Define usage line and long description.
    addUsageLine(parser, "-nm \"\\fINAMES.dmp\\fP\" -nd \"\\fINODES.dmp\\fP\" -o \"\\fISLIMM.sldb\\fP\" [\\fIOPTIONS\\fP] \"\\fIFASTA\\fP\" \"\\fIACCESSION2TAXAID\\fP\"  [\\fIACCESSION2TAXAID_2 ...\\fP]");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "FASTA FILE"));
    setValidValues(parser, 0, SeqFileIn::getFileExtensions());
    setHelpText(parser, 0, "A multi-fasta file used as a reference for mapping");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "ACCESSION2TAXAID MAP FILES", true));
    setHelpText(parser, 1, "one ore more accession to taxa id mapping files dowloaded from ncbi (separated by space.)");

    // The output file argument.
    addOption(parser, ArgParseOption("o", "output-file", "The path to the output file (default slimm_db.sldb)",
                                     ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "output-file", ".sldb");
    setDefaultValue(parser, "output-file", options.output_path);

    addOption(parser, ArgParseOption("nm", "names", "NCBI's names.dmp file which contains the mapping of taxaid to name",
                             ArgParseArgument::INPUT_FILE));
    setRequired(parser, "names");

    addOption(parser, ArgParseOption("nd", "nodes", "NCBI's nodes.dmp file which contains the taxonomic tree.",
                             ArgParseArgument::INPUT_FILE));
    setRequired(parser, "nodes");

    addOption(parser, ArgParseOption("b", "batch", "maximum number of mapping to load to memory. (default=1000000)",
                             ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "batch", options.batch);

    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
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

    getArgumentValue(options.fasta_path, parser, 0);
    uint32_t acc__taxaid_count = getArgumentValueCount(parser, 1);
    options.ac__taxid_paths.resize(acc__taxaid_count);

    for (uint32_t i = 0; i < acc__taxaid_count; ++i)
        getArgumentValue(options.ac__taxid_paths[i], parser, 1, i);    

    getOptionValue(options.names_path, parser, "names");
    getOptionValue(options.nodes_path, parser, "nodes");

    if (isSet(parser, "output-file"))
        getOptionValue(options.output_path, parser, "output-file");
    if (isSet(parser, "batch"))
        getOptionValue(options.batch, parser, "batch");
    if (isSet(parser, "verbose"))
        getOptionValue(options.verbose, parser, "verbose");

    return ArgumentParser::PARSE_OK;
}


// --------------------------------------------------------------------------
// Function get_accession_numbers()
// --------------------------------------------------------------------------
inline void get_accession_numbers(std::set<std::string> & accessions, arg_options const & options)
{
    std::cerr <<"[MSG] getting accessions numbers from fasta file ...\n";
    CharString id;
    IupacString seq;

    SeqFileIn fasta_file;
    if (!open(fasta_file, toCString(options.fasta_path)))
    {
        CharString msg = "Unable to open contigs File: ";
        append (msg, options.fasta_path);
        throw toCString(msg);
    }
    while(!atEnd(fasta_file))
    {
        readRecord(id, seq, fasta_file);
        accessions.insert(get_accession_id(id));
    }
    close(fasta_file);
}

// --------------------------------------------------------------------------
// Function get_batch_mappings_ac__taxid()
// --------------------------------------------------------------------------
inline bool get_batch_mappings_ac__taxid(std::unordered_map<std::string, uint32_t> & ac__taxid_map,
                                         std::ifstream & ac__taxid_stream,
                                         uint32_t const batch_size)
{
    ac__taxid_map.clear();
    uint32_t taxid = 0, lines_count = 0;
    std::string ac, line, ignore;

    while(std::getline(ac__taxid_stream, line))
    {
        std::stringstream   linestream(line);
        std::getline(linestream, ac, '\t'); // first column is accesion
        std::getline(linestream, ignore, '\t'); // skip the second column (accesion with version)
        linestream >> taxid;// third column is taxid
        ac__taxid_map[ac]=taxid;
        lines_count++;
        if (lines_count >= batch_size)
            break;
    }
    return (lines_count!=0);
}
// --------------------------------------------------------------------------
// Function get_taxid_from_accession()
// --------------------------------------------------------------------------
inline void get_taxid_from_accession(slimm_database & slimm_db,
                                     std::set<std::string> & accessions,
                                     arg_options const & options)
{
    std::cerr <<"[MSG] mapping accessions to taxaid ...\n";
    uint32_t accessions_count = accessions.size();
    uint32_t map_file_number  = 1;
    // iterate over multiple files
    for(std::string map_path : options.ac__taxid_paths)
    {
        if (accessions.size() == 0) // if all accesions are accounted for
            return;
        std::unordered_map<std::string, uint32_t> ac__taxid_map;
        std::ifstream ac__taxid_stream(map_path);

        // iterate over a batch of mappings: for memory sake
        uint32_t iter_number  = 1;
        while(get_batch_mappings_ac__taxid(ac__taxid_map, ac__taxid_stream, options.batch))
        {
            if (accessions.size() == 0) // if all accesions are accounted for
                return;
            if (options.verbose)
            {
                std::cerr << "[VERBOSE MSG] mapping file: ["<< map_file_number <<"/"<< options.ac__taxid_paths.size() << "]\t";
                std::cerr << "iter: [" << iter_number << "]\t";
                std::cerr << "accessions left: ["<< accessions.size() << "/" << accessions_count <<"]\n";
                ++iter_number;
            }
            // iterate over the remaining accessions to get
            for(auto ac_it=accessions.begin(); ac_it != accessions.end();)
            {
                auto ac_pos = ac__taxid_map.find(*ac_it);
                if(ac_pos != ac__taxid_map.end())
                {
                    //insert the found accessions in to the DB
                    slimm_db.ac__taxid[*ac_it] = std::vector<uint32_t>(LINAGE_LENGTH, 0);
                    slimm_db.ac__taxid[*ac_it][0] = ac_pos->second;

                    //remove found accessions form the set
                    ac_it = accessions.erase(ac_it);
                }
                else
                {
                    ++ac_it;
                }
            }
        }
        ac__taxid_stream.close();
        ++map_file_number;
    }

    // some accessions are still not mapped
    std::cerr <<"[WARNING!] The following accessions were not mapped to taxaid.\n";
    for(auto ac_it=accessions.begin(); ac_it != accessions.end(); ++ac_it)
    {
        std::cerr << *ac_it << ", ";
    }
    std::cerr <<"\nTry including the dead ACCESSION2TAXAID MAP FILE (e.g. dead_nucl.accession2taxid)\n";
}

// --------------------------------------------------------------------------
// Function fill_name_taxid_linage()
// --------------------------------------------------------------------------
inline void fill_name_taxid_linage(slimm_database & slimm_db, arg_options const & options)
{
    std::cerr <<"[MSG] loading nodes and names mappings from files ...\n";
    std::unordered_map<uint32_t, std::tuple<taxa_ranks, uint32_t> > taxid__parent;
    std::unordered_map<uint32_t, std::string>                       taxid__name;

    std::ifstream taxid__parent_stream(options.nodes_path);
    std::ifstream taxid__name_stream(options.names_path);

    uint32_t taxid=0, parent_taxid = 0;
    std::string line, rank, ignore, name;

    while(std::getline(taxid__parent_stream, line))
    {
        std::stringstream   linestream(line);
        linestream >> taxid;// first column is taxid
        std::getline(linestream, ignore, '\t'); //skip |
        std::getline(linestream, ignore, '\t'); //skip |
        linestream >> parent_taxid; // third column is parent_taxid
        std::getline(linestream, ignore, '\t'); //skip |
        std::getline(linestream, ignore, '\t'); //skip |
        std::getline(linestream, rank, '\t'); //fifth column is rank
        taxid__parent[taxid] = std::make_tuple(to_taxa_ranks(rank), parent_taxid);
    }

    taxid__parent_stream.close();

    while(std::getline(taxid__name_stream, line))
    {
        auto pos = line.find("scientific name", 0);
        if(pos != std::string::npos)
        {
            std::stringstream   linestream(line);
            linestream >> taxid;// first column is taxid
            std::getline(linestream, ignore, '\t'); //skip |
            std::getline(linestream, ignore, '\t'); //skip |
            std::getline(linestream, name, '\t'); // 2nd column is name
            taxid__name[taxid] = name;
        }
    }
    taxid__name_stream.close();

    std::cerr <<"[MSG] getting taxonomic linages and resolving names ...\n";
    for(auto ac__taxid_it=slimm_db.ac__taxid.begin(); ac__taxid_it != slimm_db.ac__taxid.end(); ++ac__taxid_it)
    {
        uint32_t tid = ac__taxid_it->second[0];
        slimm_db.taxid__name[tid] = std::make_tuple(strain_lv, taxid__name[tid]);

        while (tid != 1)
        {
            auto tid_pos = taxid__parent.find(tid);
            if(tid_pos == taxid__parent.end())
                break;

            taxa_ranks current_rank = std::get<0>(tid_pos->second);
            if (current_rank >= species_lv && current_rank <= superkingdom_lv)
            {
                ac__taxid_it->second[current_rank] = tid;
                slimm_db.taxid__name[tid] = std::make_tuple(current_rank, taxid__name[tid]);
            }
            tid = std::get<1>(tid_pos->second);
        }
    }
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

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // get the accession numbers from the fasta file
    std::set<std::string> accessions;
    get_accession_numbers(accessions, options);

    slimm_database slimm_db;
    // get the taxid from accession numbers
    get_taxid_from_accession(slimm_db, accessions, options);
    fill_name_taxid_linage(slimm_db, options);
    save_slimm_database(slimm_db, options.output_path);

//
//    std::vector<uint32_t> tids = slimm_db.ac__taxid["NZ_CP009257.1"];
//    std::cout << "ACC: NC_004578.1 \t taxa id: " << tids[0] << "\n";
//    for (uint32_t i=0; i<tids.size(); ++i)
//    {
//        std::string r = from_taxa_ranks(static_cast<taxa_ranks>(i));
//        std::cout << r << "\t" << from_taxa_ranks(std::get<0>(slimm_db.taxid__name[tids[i]])) << "\t" << tids[i] << "\t" << std::get<1>(slimm_db.taxid__name[tids[i]])  << "\n";
//    }

    return 0;
}
