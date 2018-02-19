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

#ifndef SLIMM_H
#define SLIMM_H


using namespace seqan;

inline uint32_t findLCATaxaID(std::set<uint32_t> const & taxon_ids, TNodes const & nodes);
// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class arg_options
// ----------------------------------------------------------------------------
struct arg_options
{
    typedef std::vector<std::string>            TList;

    TList rankList = {"all",
                      "species",
                      "genus",
                      "family",
                      "order",
                      "class",
                      "phylum",
                      "superkingdom"};

    float               cov_cut_off;
    uint32_t            bin_width;
    uint32_t            min_reads;
    bool                verbose;
    bool                is_directory;
    bool                raw_output;
    std::string         rank;
    std::string         input_path;
    std::string         output_prefix;
    std::string         mapping_dir;

    arg_options() : cov_cut_off(0.99),
                    bin_width(0),
                    min_reads(0),
                    verbose(false),
                    is_directory(false),
                    raw_output(false),
                    rank("all"),
                    input_path(""),
                    output_prefix(""),
                    mapping_dir("taxonomy/") {}
};

// ----------------------------------------------------------------------------
// Class slimm
// ----------------------------------------------------------------------------
class slimm
{
public:
    //constructor with argument options
    slimm(arg_options op): options(op)
    {
        collect_bam_files();
        load_taxonomic_info();
    }

    arg_options                                         options;

    uint32_t                    current_file_index        = 0;
    uint32_t                    number_of_files           = 0;
    uint16_t                    avg_read_length           = 0;
    uint32_t                    matched_ref_length        = 0;
    uint32_t                    reference_count           = 0;
    uint32_t                    failed_by_min_read        = 0;
    uint32_t                    failed_by_min_uniq_read   = 0;
    uint32_t                    failed_byCov              = 0;
    uint32_t                    failed_byUniqCov          = 0;
    uint32_t                    hits_count                = 0;
    uint32_t                    uniq_hits_count           = 0;
    uint32_t                    matches_count             = 0;
    uint32_t                    uniq_matches_count        = 0;
    uint32_t                    uniq_matches_count2       = 0;



    TNodes                                              nodes;
    std::set<uint32_t>                                  valid_ref_taxon_ids;
    std::vector<uint32_t>                               matched_taxa;
    std::vector<reference_contig>                       references;
    std::unordered_map<uint32_t, float>                 taxon_id__abundance;
    std::unordered_map <uint32_t, std::string>          taxon_id__name;
    std::unordered_map<std::string, read_stat>          reads;
    std::unordered_map<uint32_t, uint32_t>              taxon_id__readCount;
    std::unordered_map<uint32_t, std::set<uint32_t> >   taxon_id__children;

    inline std::string current_bam_file_path()
    {
        return _input_paths[current_file_index];
    }

    inline void     analyze_alignments(BamFileIn & bam_file);
    inline float    coverage_cut_off();
    inline float    expected_coverage() const;
    inline void     filter_alignments();
    inline void     get_profiles();
    inline void     get_reads_lca_count();
    inline uint32_t min_reads();
    inline uint32_t min_uniq_reads();
    inline void     print_filter_stat();
    inline void     print_matches_stat();
    inline float    uniq_coverage_cut_off();
    inline void     write_raw_output_file();
    inline void     write_output_files();
    inline void     write_output_files(uint32_t const & read_count_at_rank, std::string const & rank);
    inline void     write_to_file(std::string & filePath, std::vector<reference_contig> & refList);

private:

    float                       _coverage_cut_off       = 0.0;
    float                       _uniq_coverage_cut_off  = 0.0;
    int32_t                     _min_uniq_reads         = -1;
    int32_t                     _min_reads              = -1;
    std::vector<std::string>    _input_paths;

    // member functions
    inline void collect_bam_files();
    inline void load_taxonomic_info();
};


inline void slimm::analyze_alignments(BamFileIn & bam_file)
{
    uint32_t queryLen = 0;
    BamAlignmentRecord record;
    while (!atEnd(bam_file))
    {
        readRecord(record, bam_file);
        if (hasFlagUnmapped(record) || record.rID == BamAlignmentRecord::INVALID_REFID)
            continue;  // Skip these records.

        uint32_t newQueryLen = length(record.seq);
        queryLen = (newQueryLen == 0) ? queryLen : newQueryLen;
        uint32_t center_position =  std::min(record.beginPos + (queryLen/2),references[record.rID].length);
        uint32_t relativeBinNo = center_position/options.bin_width;

        // maintain read properties under slimm.reads
        std::string readName = toCString(record.qName);
        if(hasFlagFirst(record))
            append(readName, ".1");
        else if(hasFlagLast(record))
            append(readName, ".2");
        
        // if there is no read with readName this will create one.
        reads[readName].add_target(record.rID, relativeBinNo);
        reads[readName].len = queryLen;
        ++hits_count;
    }

    if (hits_count != 0)
    {
        uint32_t concatQLength = 0;
        for (auto it= reads.begin(); it != reads.end(); ++it)
        {
            concatQLength += it->second.len;
            if(it->second.is_uniq(matched_taxa))
            {
                uint32_t reference_id = it->second.targets[0].reference_id;
                it->second.refs_length_sum += references[reference_id].length;
                ++uniq_matches_count;

                size_t pos_count = (it->second.targets[0]).positions.size();
                references[reference_id].reads_count += pos_count;
                it->second.refs_length_sum += references[reference_id].length;
                for (size_t j=0; j < pos_count; ++j)
                {
                    uint32_t bin_number = (it->second.targets[0]).positions[j];
                    ++references[reference_id].cov.bins_height[bin_number];
                }
                references[reference_id].uniq_reads_count += 1;
                uniq_hits_count += 1;
                ++references[reference_id].uniq_cov.bins_height[(it->second.targets[0]).positions[0]];
            }
            else
            {
                size_t len = it->second.targets.size();
                for (size_t i=0; i < len; ++i)
                {

                    uint32_t reference_id = it->second.targets[i].reference_id;
                    it->second.refs_length_sum += references[reference_id].length;

                    // ***** all of the matches in multiple pos will be counted *****
                    references[reference_id].reads_count += (it->second.targets[i]).positions.size();
                    for (auto bin_number : (it->second.targets[i]).positions)
                    {
                        ++references[reference_id].cov.bins_height[bin_number];
                    }
                }
            }
        }
        matches_count = reads.size();

        avg_read_length = concatQLength/matches_count;
        float totalAb = 0.0;
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].reads_count > 0)
            {
                ++reference_count;
                matched_ref_length += references[i].length;
                references[i].relative_abundance = float(references[i].reads_count * 100)/hits_count;
                totalAb += references[i].relative_abundance/references[i].length;
            }
            else
            {
                references[i].relative_abundance = 0.0;
            }
        }
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].reads_count > 0)
            {
                references[i].relative_abundance = (references[i].relative_abundance * 100) / (totalAb*references[i].length);
            }
        }

        totalAb = 0.0;
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].uniq_reads_count > 0)
            {
                references[i].uniq_relative_abundance = float(references[i].uniq_reads_count * 100)/uniq_hits_count;
                totalAb += references[i].uniq_relative_abundance/references[i].length;
            }
            else
            {
                references[i].uniq_relative_abundance = 0.0;
            }
        }
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].uniq_reads_count > 0)
            {
                references[i].uniq_relative_abundance = (references[i].uniq_relative_abundance * 100) / (totalAb*references[i].length);
            }

        }
    }
}

//collect the sam files to process
inline void slimm::collect_bam_files()
{
    number_of_files = 1;
    if (options.is_directory)
    {
        _input_paths = get_bam_files_in_directory(options.input_path);
        number_of_files = length(_input_paths);
        if (options.verbose)
            std::cerr << number_of_files << " SAM/BAM Files found under the directory: " << options.input_path << "!\n";
    }
    else
    {
        if (is_file(toCString(options.input_path)))
            _input_paths.push_back(options.input_path);
        else
        {
            std::cerr << options.input_path << " is not a file use -d option for a directory.\n";
            exit(1);
        }
    }
}

float slimm::coverage_cut_off()
{
    if (_coverage_cut_off == 0.0 && options.cov_cut_off < 1.0)
    {
        std::vector<float> covs = {};
        covs.reserve(length(references));
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].uniq_reads_count > 0)
            {
                covs.push_back(references[i].cov_percent());
            }
        }
        _coverage_cut_off = get_quantile_cut_off<float>(covs, options.cov_cut_off);
    }
    return _coverage_cut_off;
}

float slimm::expected_coverage() const
{
    return float(avg_read_length * matches_count) / matched_ref_length;
}

inline void slimm::filter_alignments()
{
    uint32_t reference_count = length(references);
    for (uint32_t i=0; i<reference_count; ++i)
    {

        if (references[i].reads_count == 0)
            continue;
        if (
            references[i].cov_percent() >= coverage_cut_off() &&
            references[i].uniq_cov_percent() >= uniq_coverage_cut_off() &&
            true
            )
            valid_ref_taxon_ids.insert(matched_taxa[i]);
        else
        {
            if(references[i].uniq_cov_percent() < uniq_coverage_cut_off())
            {
                ++failed_byUniqCov;
            }
            if(references[i].reads_count < options.min_reads)
            {
                ++failed_by_min_read;
            }
            if(references[i].cov_percent() < coverage_cut_off())
            {
                ++failed_byCov;
            }
        }
    }

    for (auto it= reads.begin(); it != reads.end(); ++it)
    {
        it->second.update(matched_taxa, valid_ref_taxon_ids, references);
        if(it->second.is_uniq(matched_taxa))
        {
            uint32_t reference_id = (it->second.targets[0]).reference_id;
            references[reference_id].uniq_reads_count2 += 1;
            uniq_matches_count2 += 1;
            uint32_t bin_number = (it->second.targets[0]).positions[0];
            ++references[reference_id].uniq_cov2.bins_height[bin_number];
        }
    }
}

// get taxonomic profiles from the sam/bam 
inline void slimm::get_profiles()
{
    Timer<>  stop_watch;

    BamFileIn bam_file;
    BamHeader bam_header;

    std::cerr   << "\nReading " << current_file_index + 1 << " of " << number_of_files << " files ... ("
                << get_file_name(current_bam_file_path()) << ")\n"
                <<"=================================================================\n";

    if (read_bam_file(bam_file, bam_header, current_bam_file_path()))
    {
        CharString first_ref_name = contigNames(context(bam_file))[0];
        
        // Determine taxa id position
        uint32_t taxon_id_pos = get_taxon_id_pos(first_ref_name);

        //if bin_width is not given use avg read length
        if (options.bin_width == 0) 
            options.bin_width = get_bin_width_from_sample(bam_file);

        //reset the bam_file to the first recored by closing and reopening
        close(bam_file);

        read_bam_file(bam_file, bam_header, current_bam_file_path());

        StringSet<CharString>   reference_names = contigNames(context(bam_file));
        StringSet<uint32_t>     refLengths;
        refLengths = contigLengths(context(bam_file));

        uint32_t references_count = length(reference_names);
        references.resize(references_count);

        std::cerr<<"Intializing coverages for all reference genome ... ";
        // Intialize coverages for all genomes

        for (uint32_t i=0; i < references_count; ++i)
        {
            reference_contig current_ref(reference_names[i], refLengths[i], options.bin_width, taxon_id_pos);
            references[i] = current_ref;
            matched_taxa.push_back(current_ref.taxon_id);
        }
        std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;


        std::cerr<<"Analysing alignments, reads and references ....... ";
        analyze_alignments(bam_file);
        std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;
        if (hits_count > 0)
        {
            // Set the minimum reads to 10k-th of the total number of matched reads if not set by the user
            if (options.min_reads == 0)
              options.min_reads = 1 + ((matches_count - 1) / 10000);
            if (options.verbose)
                print_matches_stat();

            std::cerr   << "Filtering unlikely sequences ..................... ";
            filter_alignments();
            std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;

            if (options.verbose)
                print_filter_stat();

            if (options.raw_output)
            {
                std::cerr<<"Writing features to a file ....................... ";
                write_raw_output_file();
                std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;
            }

            std::cerr<<"Assigning reads to Least Common Ancestor (LCA) ... ";
            get_reads_lca_count();
            std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;

            std::cerr<<"Writing taxnomic profile(s) ...................... ";
            write_output_files();
            if (options.verbose)
                std::cerr<<"\n.................................................. ";
            std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;


        }

        else
        {
            std::cerr << "[WARNING] No mapped reads found in BAM file!" << std::endl;
        }
        
        std::cerr<<"[Done!] File took " << stop_watch.elapsed() <<" secs to process.\n";
    }
}

inline void slimm::get_reads_lca_count()
{
    // put the non-unique read to upper taxa.
    for (auto it= reads.begin(); it != reads.end(); ++it)
    {
        if(!(it->second.is_uniq(matched_taxa)))
        {
            uint32_t lcaTaxaID = 0;
            std::set<uint32_t> ref_ids = {};
            std::set<uint32_t> taxon_ids;
            size_t len = it->second.targets.size();
            for (size_t i=0; i < len; ++i)
            {
                uint32_t ref_id = (it->second.targets[i]).reference_id;
                taxon_ids.insert(matched_taxa[ref_id]);
                ref_ids.insert(ref_id);
            }
            lcaTaxaID = getLCA(taxon_ids, nodes);
            taxon_id__children[lcaTaxaID].insert(ref_ids.begin(), ref_ids.end());
            // If taxon_id already exists
            if (taxon_id__readCount.count(lcaTaxaID) == 1)
                ++taxon_id__readCount[lcaTaxaID];
            else   // first time for taxon_id
                taxon_id__readCount[lcaTaxaID] = 1;
            //add the contributing children references to the taxa
        }
    }
    //add the sum of read counts of children all ancestors of the LCA
    std::unordered_map <uint32_t, uint32_t> tID2ReadCountCopy = taxon_id__readCount;
    for (auto t2rc : tID2ReadCountCopy)
    {
        uint32_t currentTaxaID = t2rc.first;
        uint32_t readCount = tID2ReadCountCopy[currentTaxaID];
        std::set<uint32_t> ref_ids = taxon_id__children[currentTaxaID];
        while (nodes.count(currentTaxaID) == 1 && currentTaxaID != 0)
        {
            currentTaxaID = (nodes.at(currentTaxaID)).first;
            if (taxon_id__readCount.count(currentTaxaID) >= 1)
                taxon_id__readCount[currentTaxaID] += readCount;
            else
                taxon_id__readCount[currentTaxaID] = readCount;
            //add the contributing children references to the taxa
            taxon_id__children[currentTaxaID].insert(ref_ids.begin(), ref_ids.end());
        }
    }


    for (uint32_t i=0; i<length(references); ++i)
    {

        if (references[i].uniq_reads_count2 > 0)
        {
            uint32_t currentTaxaID = references[i].taxon_id;
            uint32_t uniqCount = references[i].uniq_reads_count2;
            while (nodes.count(currentTaxaID) == 1 && currentTaxaID != 0)
            {
                if (taxon_id__readCount.count(currentTaxaID) >= 1)
                    taxon_id__readCount[currentTaxaID] += uniqCount;
                else
                    taxon_id__readCount[currentTaxaID] = uniqCount;
                taxon_id__children[currentTaxaID].insert(i);
                //add the contributing children references to the taxa
                currentTaxaID = (nodes.at(currentTaxaID)).first;
            }
        }
    }
}

inline void slimm::print_filter_stat()
{
    std::cerr << "  " << length(valid_ref_taxon_ids) << " passed the threshould coverage.\n";
    std::cerr << "  " << failed_byCov << " ref's couldn't pass the coverage threshould.\n";
    std::cerr << "  " << failed_byUniqCov << " ref's couldn't pass the uniq coverage threshould.\n";
    std::cerr << "  uniquily matching reads increased from " << uniq_matches_count << " to " << uniq_matches_count2 <<"\n\n";
}

inline void slimm::print_matches_stat()
{
    std::cerr << "  "   << hits_count << " records processed." << std::endl;
    std::cerr << "    " << matches_count << " matching reads" << std::endl;
    std::cerr << "    " << uniq_matches_count << " uniquily matching reads"<< std::endl;
    std::cerr << "  references with reads = " << reference_count << std::endl;
    std::cerr << "  expected bins coverage = " << expected_coverage() <<std::endl;
    std::cerr << "  bins coverage cut-off = " << coverage_cut_off() << " (" << options.cov_cut_off <<" quantile)\n";
    std::cerr << "  uniq bins coverage cut-off = " << uniq_coverage_cut_off() << " (" << options.cov_cut_off <<" quantile)\n\n";
}

// load taxonomic information from slimmDB
inline void slimm::load_taxonomic_info()
{
    Timer<>  stop_watch;
    std::stringstream ss;
    std::ofstream sam_extract_file;

    std::string nodes_path = options.mapping_dir + "/nodes.dmp";
    std::string names_path = options.mapping_dir + "/names.dmp";

    std::cerr<<"Loading taxon_id to name mapping ................. ";
    taxon_id__name =  load_int__string_map(names_path);
    std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;

    std::cerr<<"Loading node mapping ............................. ";
    nodes = load_node_maps(nodes_path);
    std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;
}

uint32_t slimm::min_reads()
{
    if (options.cov_cut_off == 1.0)
        _min_reads = 0;
    if (_min_reads == -1)
    {
        std::vector<int> counts = {};
        counts.reserve(length(references));
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].reads_count > 0)
            {
                counts.push_back(references[i].reads_count);
            }
        }
        _min_reads = get_quantile_cut_off(counts, options.cov_cut_off);
    }
    return _min_reads;
}

uint32_t slimm::min_uniq_reads()
{
    if (options.cov_cut_off == 1.0)
        _min_uniq_reads = 0;
    if (_min_uniq_reads == -1)
    {
        std::vector<int> uniqCounts = {};
        uniqCounts.reserve(length(references));
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].uniq_reads_count > 0)
            {
                uniqCounts.push_back(references[i].uniq_reads_count);
            }
        }
        _min_uniq_reads = get_quantile_cut_off(uniqCounts, options.cov_cut_off);
    }
    return _min_uniq_reads;
}

float slimm::uniq_coverage_cut_off()
{
    if (_uniq_coverage_cut_off == 0.0 && options.cov_cut_off < 1.0)
    {
        std::vector<float> covs = {};
        covs.reserve(length(references));
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].uniq_reads_count > 0)
            {
                covs.push_back(references[i].uniq_cov_percent());
            }
        }
        _uniq_coverage_cut_off = get_quantile_cut_off<float>(covs, options.cov_cut_off);
    }
    return _uniq_coverage_cut_off;
}

inline void slimm::write_output_files(uint32_t const & read_count_at_rank,
                                      std::string const & rank)
{

    std::string tsvFile = get_tsv_file_name(toCString(options.output_prefix), current_bam_file_path(), rank);
    std::ofstream abundunceFile;
    abundunceFile.open(tsvFile);

    uint32_t unknownReads = matches_count-read_count_at_rank;
    uint32_t count = 0;
    uint32_t faild_count = 0;
    float faildAbundunce = 0.0;
    std::unordered_map <uint32_t, float> cladeCov;
    std::unordered_map <uint32_t, float> cladeAbundance;
    count = 1;
    abundunceFile<<"No.\tName\tTaxid\tNoOfReads\tRelativeAbundance\tContributers\tbins_coverage\n";
    for (auto tID : taxon_id__readCount) {
        if (rank == nodes[tID.first].second)
        {
            uint32_t cLength = 0;
            uint32_t noOfContribs = 0;
            std::set<uint32_t>::iterator it;
            for (it=taxon_id__children.at(tID.first).begin();
                 it!=taxon_id__children.at(tID.first).end(); ++it)
            {
                cLength += references[*it].length;
                ++noOfContribs;
            }
            cLength = cLength/noOfContribs;
            float cov = float(tID.second * avg_read_length)/cLength;

            float relative_abundance = float(tID.second)/(matches_count) * 100;
            std::unordered_map <uint32_t, std::string>::const_iterator it2 =
            taxon_id__name.find (tID.first);
            if (relative_abundance == 0.001 || cov < coverage_cut_off())
            {
                faildAbundunce += relative_abundance;
                ++faild_count;
                continue;
            }
            seqan::CharString candidateName = "Organism name not found";
            if (it2 != taxon_id__name.end())
                candidateName = (taxon_id__name.at(tID.first));
            abundunceFile   << count << "\t"
            << candidateName << "\t"
            << tID.first << "\t"
            << tID.second << "\t"
            << relative_abundance << "\t"
            << taxon_id__children.at(tID.first).size() << "\t"
            << cov << "\n";
            ++count;
        }

    }
    float unknownAbundance = float(unknownReads)/ matches_count * 100 + faildAbundunce;

    abundunceFile   << count << "\tunknown_"<< rank<< "(multiple)" << "\t0\t"
    << unknownReads << "\t" << unknownAbundance << "\t0\t0.0\n";
    abundunceFile.close();
    if (options.verbose)
    {
        std::cerr << "\n" << std::setw (15) << rank <<" level: "<< faild_count <<" bellow cutoff ("<< 0.001 <<")";
    }
}

inline void slimm::write_output_files()
{
    std::vector<std::string> const all_ranks= {"species", "genus", "family", "order", "class", "phylum", "superkingdom"};
    std::vector<std::string> considered_ranks = {options.rank};
    if(options.rank == "all")
    {
        considered_ranks.resize(7);
        considered_ranks = all_ranks;
    }

    // calculate the total number of reads matching uniquily at that species level.
    std::map<std::string, uint32_t>  read_count_at_rank = { {"all", 0},
                                                            {"species", 0},
                                                            {"genus", 0},
                                                            {"family", 0},
                                                            {"order", 0},
                                                            {"class", 0},
                                                            {"phylum", 0},
                                                            {"superkingdom", 0}};
    for (auto tID : taxon_id__readCount)
    {
        if (read_count_at_rank.find(nodes[tID.first].second) != read_count_at_rank.end() )
        {
            read_count_at_rank[nodes[tID.first].second] +=  tID.second ;     // found
        }
    }

    for (std::string rank : considered_ranks)
    {
        write_output_files(read_count_at_rank[rank], rank);
    }
}

inline void slimm::write_raw_output_file()
{
    std::string tsvFile = get_tsv_file_name(options.output_prefix, current_bam_file_path());
    write_to_file(tsvFile, references);
}

inline void slimm::write_to_file(std::string & filePath, std::vector<reference_contig> & refList)
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
    "none_zero_bin_count\t"
    "none_zero_bin_countUniq\t"
    "none_zero_bin_countUniq2\t"


    "bins_coverageDepth\t"

    "Uniqbins_coverageDepth\t"
    "Uniqbins_coverageDepth2\t"

    "MappingError\t"
    "bins_coveragePercentage\t"
    "Uniquebins_coveragePercentage\t"
    "Uniquebins_coveragePercentage2\n";

    uint32_t current = 0;
    uint32_t reference_count = length(refList);
    for (uint32_t i=0; i < reference_count; ++i)
    {
        current ++;
        reference_contig current_ref = refList[i];
        CharString candidateName = current_ref.reference_name;
        std::unordered_map <uint32_t, std::string>::const_iterator it = taxon_id__name.find(current_ref.taxon_id);
        if (it != taxon_id__name.end())
            candidateName = (taxon_id__name.at(current_ref.taxon_id));

        features_file   << current << "\t"
        << candidateName << "\t"
        << current_ref.taxon_id << "\t"
        << current_ref.reads_count << "\t"
        << current_ref.relative_abundance << "\t"
        << current_ref.uniq_relative_abundance << "\t"
        << current_ref.uniq_relative_abundance2 << "\t"
        << current_ref.length << "\t"
        << current_ref.uniq_reads_count << "\t"
        << current_ref.uniq_reads_count2 << "\t"

        << current_ref.cov.number_of_bins << "\t"
        << current_ref.cov.none_zero_bin_count() << "\t"
        << current_ref.uniq_cov.none_zero_bin_count() << "\t"
        << current_ref.uniq_cov2.none_zero_bin_count() << "\t"


        << current_ref.cov_depth() << "\t"

        << current_ref.uniq_cov_depth() << "\t"
        << current_ref.uniq_cov_depth2() << "\t"

        << "NA"<< "\t"
        << current_ref.cov_percent() << "\t"
        << current_ref.uniq_cov_percent() << "\t"
        << current_ref.uniq_cov_percent2() << "\n";
    }
    features_file.close();
}




inline int get_taxonomic_profile(arg_options & options)
{
    // slimm object
    Timer<>  stop_watch;
    uint32_t total_hits_count = 0;
    slimm slimm1(options);
    for (uint32_t n=0; n < slimm1.number_of_files; ++n)
    {
        slimm1.current_file_index = n;
        slimm1.get_profiles();
        total_hits_count += slimm1.hits_count;
    }

    std::string output_directory = get_directory(options.output_prefix);

    std::cerr << "\n*****************************************************************\n";
    std::cerr << total_hits_count << " SAM/BAM alignment records are proccessed.\n";
    std::cerr << "Taxonomic profiles are written to: \n   " << output_directory <<"\n";
    std::cerr << "Total time elapsed: " << stop_watch.elapsed() <<" secs\n";

    return 0;
}


#endif /* SLIMM_H */
