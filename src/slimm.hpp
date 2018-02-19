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
    std::string         database_path;

    arg_options() : cov_cut_off(0.99),
                    bin_width(0),
                    min_reads(0),
                    verbose(false),
                    is_directory(false),
                    raw_output(false),
                    rank("all"),
                    input_path(""),
                    output_prefix(""),
                    database_path("") {}
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
        get_considered_ranks();
        load_slimm_database(db, options.database_path);
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



    slimm_database                                      db;
    std::set<uint32_t>                                  valid_ref_taxon_ids;
    std::vector<uint32_t>                               matched_taxa;
    std::vector<taxa_ranks>                             considered_ranks;
    std::vector<reference_contig>                       references;
    std::unordered_map<uint32_t, float>                 taxon_id__abundance;
//    std::unordered_map <uint32_t, std::string>          taxon_id__name;
    std::unordered_map<std::string, read_stat>          reads;
    std::unordered_map<uint32_t, uint32_t>              taxon_id__read_count;
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
    inline void     write_output_files(uint32_t const & read_count_at_rank, taxa_ranks const & rank);
    inline void     write_to_file(std::string & filePath, std::vector<reference_contig> & refList);

private:

    float                       _coverage_cut_off       = 0.0;
    float                       _uniq_coverage_cut_off  = 0.0;
    int32_t                     _min_uniq_reads         = -1;
    int32_t                     _min_reads              = -1;
    std::vector<std::string>    _input_paths;

    // member functions
    inline void collect_bam_files();
    inline void get_considered_ranks();
    inline void load_taxonomic_info();
};


inline void slimm::analyze_alignments(BamFileIn & bam_file)
{
    uint32_t query_len = 0;
    BamAlignmentRecord record;
    while (!atEnd(bam_file))
    {
        readRecord(record, bam_file);
        if (hasFlagUnmapped(record) || record.rID == BamAlignmentRecord::INVALID_REFID)
            continue;  // Skip these records.

        uint32_t new_query_len = length(record.seq);
        query_len = (new_query_len == 0) ? query_len : new_query_len;
        uint32_t center_position =  std::min(record.beginPos + (query_len/2), references[record.rID].length);
        uint32_t relative_bin_no = center_position/options.bin_width;

        // maintain read properties under slimm.reads
        std::string read_name = toCString(record.qName);
        if(hasFlagFirst(record))
            append(read_name, ".1");
        else if(hasFlagLast(record))
            append(read_name, ".2");
        
        // if there is no read with read_name this will create one.
        reads[read_name].add_target(record.rID, relative_bin_no);
        reads[read_name].len = query_len;
        ++hits_count;
    }


    if (hits_count == 0)
        return;

    uint32_t query_length_sum = 0;
    for (auto it= reads.begin(); it != reads.end(); ++it)
    {
        query_length_sum += it->second.len;
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

    avg_read_length = query_length_sum/matches_count;
    float totalAb = 0.0;
    for (uint32_t i=0; i<length(references); ++i)
    {
        if (references[i].reads_count > 0)
        {
            ++reference_count;
            matched_ref_length += references[i].length;
            references[i].abundance = float(references[i].reads_count * 100)/hits_count;
            totalAb += references[i].abundance/references[i].length;
        }
        else
        {
            references[i].abundance = 0.0;
        }
    }
    for (uint32_t i=0; i<length(references); ++i)
    {
        if (references[i].reads_count > 0)
        {
            references[i].abundance = (references[i].abundance * 100) / (totalAb*references[i].length);
        }
    }

    totalAb = 0.0;
    for (uint32_t i=0; i<length(references); ++i)
    {
        if (references[i].uniq_reads_count > 0)
        {
            references[i].uniq_abundance = float(references[i].uniq_reads_count * 100)/uniq_hits_count;
            totalAb += references[i].uniq_abundance/references[i].length;
        }
        else
        {
            references[i].uniq_abundance = 0.0;
        }
    }
    for (uint32_t i=0; i<length(references); ++i)
    {
        if (references[i].uniq_reads_count > 0)
        {
            references[i].uniq_abundance = (references[i].uniq_abundance * 100) / (totalAb*references[i].length);
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
    for (uint32_t i=0; i < reference_count; ++i)
    {
        if (references[i].reads_count == 0)
            continue;
        if (references[i].cov_percent() >= coverage_cut_off() &&
            references[i].uniq_cov_percent() >= uniq_coverage_cut_off() && true)
        {
            valid_ref_taxon_ids.insert(db.ac__taxid[references[i].accession][0]);
        }
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
        //if bin_width is not given use avg read length
        if (options.bin_width == 0) 
            options.bin_width = get_bin_width_from_sample(bam_file);

        //reset the bam_file to the first recored by closing and reopening
        close(bam_file);

        read_bam_file(bam_file, bam_header, current_bam_file_path());

        StringSet<CharString>    contig_names = contigNames(context(bam_file));
        StringSet<uint32_t>      refLengths;
        refLengths = contigLengths(context(bam_file));

        uint32_t references_count = length(contig_names);
        references.resize(references_count);

        std::cerr<<"Intializing coverages for all reference genome ... ";
        // Intialize coverages for all genomes

        for (uint32_t i=0; i < references_count; ++i)
        {
            std::string accession(toCString(contig_names[i]));
            reference_contig current_ref(accession, db.ac__taxid[accession][0], refLengths[i], options.bin_width);
            references[i] = current_ref;
            matched_taxa.push_back(db.ac__taxid[accession][0]);
        }
        std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;

        std::cerr<<"Analysing alignments, reads and references ....... ";
        analyze_alignments(bam_file);
        std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;
        if (hits_count == 0)
        {
            std::cerr << "[WARNING] No mapped reads found in BAM file!" << std::endl;
            return;
        }

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

        std::cerr<<"[Done!] File took " << stop_watch.elapsed() <<" secs to process.\n";
    }
}

inline void slimm::get_considered_ranks()
{
    if(options.rank == "all")
    {
        for(uint32_t i=1; i<8; ++i)
            considered_ranks.push_back(static_cast<taxa_ranks>(i));
    }
    else
    {
        considered_ranks.push_back(to_taxa_ranks(options.rank));
    }
}
inline void slimm::get_reads_lca_count()
{
    // put the non-unique read to upper taxa.
    for (auto it= reads.begin(); it != reads.end(); ++it)
    {
        if(!(it->second.is_uniq(matched_taxa)))
        {
            uint32_t lca_taxa_id = 0;
            std::set<uint32_t> ref_ids = {};
            std::set<uint32_t> taxon_ids;
            size_t len = it->second.targets.size();
            for (size_t i=0; i < len; ++i)
            {
                uint32_t ref_id = (it->second.targets[i]).reference_id;
                taxon_ids.insert(matched_taxa[ref_id]);
                ref_ids.insert(ref_id);
            }
            lca_taxa_id = get_lca(taxon_ids, db);
            taxon_id__children[lca_taxa_id].insert(ref_ids.begin(), ref_ids.end());
            // If taxon_id already exists
            if (taxon_id__read_count.count(lca_taxa_id) == 1)
                ++taxon_id__read_count[lca_taxa_id];
            else   // first time for taxon_id
                taxon_id__read_count[lca_taxa_id] = 1;
            //add the contributing children references to the taxa
        }
    }


    //add the sum of read counts of children to all ancestors of the LCA
    std::unordered_map <uint32_t, uint32_t> taxon_id__read_count_cp = taxon_id__read_count;
    uint32_t contributer_taxa_id = 0,  reciever_taxa_id = 0, read_count = 0, uniqCount = 0;
    for(auto ac__taxid_it=db.ac__taxid.begin(); ac__taxid_it != db.ac__taxid.end(); ++ac__taxid_it)
    {
        for (uint32_t i=0; i<LINAGE_LENGTH-1; ++i)
        {
            contributer_taxa_id = ac__taxid_it->second[i];
            if (taxon_id__read_count_cp.find(contributer_taxa_id) != taxon_id__read_count_cp.end())
            {
                std::set<uint32_t> ref_ids = taxon_id__children[contributer_taxa_id];
                read_count = taxon_id__read_count_cp[contributer_taxa_id];
                for (uint32_t j=i+1; i<LINAGE_LENGTH; ++i)
                {
                    reciever_taxa_id = ac__taxid_it->second[j];
                    if (taxon_id__read_count.count(reciever_taxa_id) >= 1)
                        taxon_id__read_count[reciever_taxa_id] += read_count;
                    else
                        taxon_id__read_count[reciever_taxa_id] = read_count;
                    //add the contributing children references to the taxa
                    taxon_id__children[reciever_taxa_id].insert(ref_ids.begin(), ref_ids.end());
                }

            }

        }
    }


    for (uint32_t i=0; i<length(references); ++i)
    {
        if (references[i].uniq_reads_count2 > 0)
        {
            std::vector<uint32_t> linage = db.ac__taxid[references[i].accession];
            uniqCount = references[i].uniq_reads_count2;
            std::set<uint32_t> ref_ids = taxon_id__children[linage[0]];
            for (uint32_t j=1; j<LINAGE_LENGTH; ++j)
            {
                reciever_taxa_id = linage[j];
                if (taxon_id__read_count.count(reciever_taxa_id) >= 1)
                    taxon_id__read_count[reciever_taxa_id] += uniqCount;
                else
                    taxon_id__read_count[reciever_taxa_id] = uniqCount;
                //add the contributing children references to the taxa
                taxon_id__children[reciever_taxa_id].insert(i);
                taxon_id__children[reciever_taxa_id].insert(ref_ids.begin(), ref_ids.end());
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
//inline void slimm::load_taxonomic_info()
//{
//    Timer<>  stop_watch;
//    std::stringstream ss;
//    std::ofstream sam_extract_file;
//
//    std::string nodes_path = options.mapping_dir + "/nodes.dmp";
//    std::string names_path = options.mapping_dir + "/names.dmp";
//
//    std::cerr<<"Loading taxon_id to name mapping ................. ";
//    taxon_id__name =  load_int__string_map(names_path);
//    std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;
//
//    std::cerr<<"Loading node mapping ............................. ";
//    nodes = load_node_maps(nodes_path);
//    std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;
//}

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

inline void slimm::write_output_files(uint32_t const & read_count_at_rank, taxa_ranks  const & rank)
{

    std::string abundunce_tsv_path = get_tsv_file_name(toCString(options.output_prefix), current_bam_file_path(), from_taxa_ranks(rank));
    std::ofstream abundunce_stream;
    abundunce_stream.open(abundunce_tsv_path);

    uint32_t unknown_reads_count = matches_count - read_count_at_rank;
    uint32_t count = 0;
    uint32_t faild_count = 0;
    float faild_abundunce = 0.0;
    std::unordered_map <uint32_t, float> clade_cov;
    std::unordered_map <uint32_t, float> clade_abundance;
    count = 1;

    abundunce_stream << "taxa_level\ttaxa_id\tlinage\tname\tabundance\tread_count\tchildren_count\n";
    for (auto t_id : taxon_id__read_count)
    {
        if (rank == std::get<0>(db.taxid__name[t_id.first]))
        {
            std::set<uint32_t>::iterator it;
            uint32_t genome_Length = 0;
            uint32_t children_count = 0;
            std::string last_child_acc = "";
            for (auto child : taxon_id__children.at(t_id.first))
            {
                genome_Length += references[child].length;
                last_child_acc = references[child].accession;
                ++children_count;
            }
            genome_Length = genome_Length/children_count;
            float cov = float(t_id.second * avg_read_length)/genome_Length;

            float abundance = float(t_id.second)/(matches_count) * 100;

            if (abundance == 0.001 || cov < coverage_cut_off())
            {
                faild_abundunce += abundance;
                ++faild_count;
                continue;
            }

            std::vector<uint32_t> linage = db.ac__taxid[last_child_acc];
            std::string linage_str = "";
            for (uint32_t i=LINAGE_LENGTH-1; i > rank; --i)
            {
                linage_str += " >";
                linage_str += std::get<1>(db.taxid__name[linage[i]]);
            }

            std::string candidate_name = std::get<1>(db.taxid__name[t_id.first]);
            if (candidate_name == "")
                candidate_name = "no_name_" + from_taxa_ranks(rank);
            abundunce_stream << from_taxa_ranks(rank) << "\t" << t_id.first << "\t" << linage_str << "\t" << candidate_name << "\t";
            abundunce_stream << abundance << "\t" << t_id.second << "\t" << children_count << "\n";
            ++count;
        }

    }
    float unknownAbundance = float(unknown_reads_count)/ matches_count * 100 + faild_abundunce;

    abundunce_stream << from_taxa_ranks(rank) << "\t" << "0" << "\t" << "-" << "\t" << "\tunknown_"<< from_taxa_ranks(rank) << " (multiple)" << "\t";
    abundunce_stream << unknownAbundance << "\t" << unknown_reads_count << "\t" << "0" << "\n";
    abundunce_stream.close();
    if (options.verbose)
    {
        std::cerr << "\n" << std::setw (15) << rank <<" level: "<< faild_count <<" bellow cutoff ("<< 0.001 <<")";
    }
}

inline void slimm::write_output_files()
{
    // calculate the total number of reads matching uniquily at that species level.
    std::map<taxa_ranks, uint32_t>  read_count_at_rank = {  {strain_lv, 0},
                                                            {species_lv, 0},
                                                            {genus_lv, 0},
                                                            {family_lv, 0},
                                                            {order_lv, 0},
                                                            {class_lv, 0},
                                                            {phylum_lv, 0},
                                                            {superkingdom_lv, 0}};

    for (auto t_id : taxon_id__read_count)
    {
        if (read_count_at_rank.find(std::get<0>(db.taxid__name[t_id.first])) != read_count_at_rank.end() )
        {
            read_count_at_rank[std::get<0>(db.taxid__name[t_id.first])] +=  t_id.second ;     // found
        }
    }

    for (taxa_ranks rank : considered_ranks)
    {
        write_output_files(read_count_at_rank[rank], rank);
    }
}

inline void slimm::write_raw_output_file()
{
    std::string tsv_path = get_tsv_file_name(options.output_prefix, current_bam_file_path());
    write_to_file(tsv_path, references);
}

inline void slimm::write_to_file(std::string & filePath, std::vector<reference_contig> & refList)
{
    std::ofstream features_file;
    features_file.open(filePath);

    features_file << "No.\t"
    "accesion\t"
    "taxaid\t"
    "name\t"
    "reads_count\t"
    "abundance\t"
    "uniq1_abundance\t"
    "uniq2_abundance\t"
    "genome_length\t"
    "uniq1_reads_count\t"
    "uniq2_reads_count\t"

    "bins_count\t"
    "bins_count(>0)\t"
    "uniq1_bins_count(>0)\t"
    "uniq2_bins_count(>0)\t"

    "coverage_depth\t"
    "uniq1_coverage_depth\t"
    "uniq2_coverage_depth\t"

    "mapping_error\t"
    "coverage(%)\t"
    "uniq1_coverage(%)\t"
    "uniq2_coverage(%)\n";

    uint32_t current = 0;
    uint32_t reference_count = length(refList);
    for (uint32_t i=0; i < reference_count; ++i)
    {
        current ++;
        reference_contig current_ref = refList[i];
        uint32_t taxa_id = db.ac__taxid[current_ref.accession][0];
        std::string candidate_name = std::get<1>(db.taxid__name[taxa_id]);
        if (candidate_name == "")
            candidate_name = "no_name_found";
        features_file   << current << "\t"
        << current_ref.accession << "\t"
        << taxa_id << "\t"
        << candidate_name << "\t"
        << current_ref.reads_count << "\t"
        << current_ref.abundance << "\t"
        << current_ref.uniq_abundance << "\t"
        << current_ref.uniq_abundance2 << "\t"
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
