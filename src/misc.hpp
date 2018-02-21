using namespace seqan;

typedef std::unordered_map <uint32_t, std::pair<uint32_t, std::string> > TNodes;
uint32_t const LINAGE_LENGTH = 8;

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <utility>

#include <cereal/types/common.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/archives/binary.hpp>

// ==========================================================================
// Classes
// ==========================================================================

enum taxa_ranks
{
    strain_lv         = 0,
    species_lv        = 1,
    genus_lv          = 2,
    family_lv         = 3,
    order_lv          = 4,
    class_lv          = 5,
    phylum_lv         = 6,
    superkingdom_lv   = 7,
    intermidiate_lv   = 8
};

taxa_ranks to_taxa_ranks(const std::string &str)
{
    if       (str == "strain")       return strain_lv;
    else if  (str == "species")      return species_lv;
    else if  (str == "genus")        return genus_lv;
    else if  (str == "family")       return family_lv;
    else if  (str == "order")        return order_lv;
    else if  (str == "class")        return class_lv;
    else if  (str == "phylum")       return phylum_lv;
    else if  (str == "superkingdom") return superkingdom_lv;
    else                             return intermidiate_lv;
}


std::string from_taxa_ranks(const taxa_ranks &rnk)
{
    if       (rnk == strain_lv)       return "strain";
    else if  (rnk == species_lv)      return "species";
    else if  (rnk == genus_lv)        return "genus";
    else if  (rnk == family_lv)       return "family";
    else if  (rnk == order_lv)        return "order";
    else if  (rnk == class_lv)        return "class";
    else if  (rnk == phylum_lv)       return "phylum";
    else if  (rnk == superkingdom_lv) return "superkingdom";
    else                              return "intermidiate";
}

struct slimm_database
{
public:
    std::unordered_map<std::string, std::vector<uint32_t> >             ac__taxid;
    std::unordered_map<uint32_t, std::tuple<taxa_ranks, std::string> >  taxid__name;

    template <class Archive>
    void save( Archive & ar ) const
    {
        ar(ac__taxid);
        ar(taxid__name);
    }

    template <class Archive>
    void load( Archive & ar )
    {
        ar(ac__taxid);
        ar(taxid__name);
    }

};

template <typename TTarget, typename TString, typename TKey = uint32_t, typename TValue = uint32_t>
TTarget load_node_maps_2(TString const & filePath)
{
    TTarget result;
    std::ifstream giMap(toCString(filePath));
    TKey key;
    TValue value;
    while (giMap >> key >> value)
    {
        result[key]=value;
    }
    giMap.close();
    return result;
}

std::unordered_map <uint32_t, std::string> load_int__string_map(std::string const & filePath)
{
    std::unordered_map <uint32_t, std::string> result;
    std::ifstream nameMap(filePath);
    uint32_t key;
    std::string value,  line;

    while(std::getline(nameMap, line))
    {
        std::stringstream   linestream(line);
        linestream >> key;
        std::getline(linestream, value, '\t');
        std::getline(linestream, value, '\t');
        result[key]=value;
    }
    nameMap.close();
    return result;
}

TNodes load_node_maps(std::string const & filePath)
{
    TNodes target;
    std::ifstream nodeMap(filePath);
    std::string   line;

    while(std::getline(nodeMap, line))
    {
        std::stringstream   linestream(line);
        uint32_t key, value1;
        std::string value2;
        linestream >> key >> value1;
        std::getline(linestream, value2, '\t');
        std::getline(linestream, value2, '\t');
        target[key]=std::pair<uint32_t, std::string> (value1, value2);
    }

    nodeMap.close();
    return target;
}

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function save_slimm_database()
// --------------------------------------------------------------------------
inline void save_slimm_database(slimm_database const & slimm_db, std::string const & output_path)
{
    std::ofstream os(output_path, std::ios::binary);
    cereal::BinaryOutputArchive out_archive( os );
    out_archive(slimm_db);
    os.close();
}

// --------------------------------------------------------------------------
// Function load_slimm_database()
// --------------------------------------------------------------------------
inline void load_slimm_database(slimm_database & slimm_db, std::string const & input_path)
{
    std::ifstream is(input_path, std::ios::binary);
    cereal::BinaryInputArchive in_archive(is);
    in_archive(slimm_db);
    is.close();
}

template <typename Type>
Type get_quantile_cut_off (std::vector<Type> v, float q)
{
    if (v.empty())
        return 0;
    Type total = std::accumulate(v.begin(), v.end(), (Type)0);
    Type cutoff = (Type)0, subTotal = (Type)0;

    std::sort (v.begin(), v.end());

    uint32_t i= v.size() - 1;
    while((float(subTotal)/total) < q && i > (Type)0)
    {
        subTotal += v[i];
        --i;
    }
    cutoff = v[i];

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
std::vector<std::string> & split(const std::string &s,
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

std::string get_file_name (const std::string& str)
{
    std::size_t found = str.find_last_of("/\\");
    return str.substr(found+1);
}

std::string get_directory (const std::string& str)
{
    std::size_t found = str.find_last_of("/\\");
    return str.substr(0,found);
}

std::string get_tsv_file_name(const std::string & output_prefix, const std::string& input_path)
{
    std::string dir_name = get_directory(output_prefix);
    std::string file_name = get_file_name(output_prefix);
    if (file_name.size() == 0)
    {
        file_name = get_file_name(input_path);

        if ((file_name.find(".sam") != std::string::npos &&
             file_name.find(".sam") == file_name.find_last_of("."))
            ||
            (file_name.find(".bam") != std::string::npos &&
             file_name.find(".bam")  == file_name.find_last_of(".")))
        {
            file_name.replace((file_name.find_last_of(".")), 4, "");
        }
    }
    return dir_name + "/" + file_name;
}

std::string get_tsv_file_name (const std::string& output_prefix, const std::string& input_path, const std::string& decor_suffix)
{
    return get_tsv_file_name(output_prefix, input_path) + decor_suffix + ".tsv";
}


// ----------------------------------------------------------------------------
// Function setDateAndVersion()
// ----------------------------------------------------------------------------

void setDateAndVersion(ArgumentParser & parser)
{
    setDate(parser, __DATE__);
    setCategory(parser, "Metagenomics");

#if defined(SEQAN_APP_VERSION)
    setVersion(parser, SEQAN_APP_VERSION);
#endif
}

// ----------------------------------------------------------------------------
// Function setDescription()
// ----------------------------------------------------------------------------

void setDescription(ArgumentParser & parser)
{
    addDescription(parser, "SLIMM  Species Level Identification of Microbes from Metagenomes");
    addDescription(parser, "See \\fI http://www.seqan.de/projects/slimm \\fP for more information.");
    addDescription(parser, "(c) Copyright 2014-2017  by Temesgen H. Dadi.");
}

typedef std::unordered_map <uint32_t, std::pair<uint32_t, std::string> > TNodes;
uint32_t getLCA(std::set<uint32_t> const & taxon_ids, std::set<uint32_t> const & valTaxaIDs, TNodes const & nodes)
{
    //consider only those under validTaxaIDs
    std::set<uint32_t> parents;
    for (auto tID : taxon_ids)
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

std::string get_accession_id(CharString const & sequence_name)
{
    typedef OrFunctor<IsWhitespace, OrFunctor<EqualsChar<'.'> , EqualsChar<'|'> > >  IsSeqNameDelim;
    StringSet <CharString> chunks;
    strSplit(chunks, sequence_name, IsSeqNameDelim());
    std::string result = toCString(chunks[0]);
    return result;
}


bool get_taxon_id(uint32_t &idPosition, CharString accession, std::string idType)
{
    StringSet <CharString> chunks;
    strSplit(chunks, accession, EqualsChar<'|'>());
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

uint32_t getLCA(std::set<uint32_t> const & taxon_ids, TNodes const & nodes)
{
    return getLCA(taxon_ids, taxon_ids, nodes);
}



uint32_t getLCA(std::vector<uint32_t> const & taxon_ids, TNodes const & nodes)
{
    std::set<uint32_t> s(taxon_ids.begin(), taxon_ids.end());
    if (s.size() == 1)
        return taxon_ids[0];
    else
        return getLCA(s, s, nodes);
}


uint32_t get_lca(std::set<uint32_t> const & taxon_ids, std::set<uint32_t> const & valid_taxon_ids, slimm_database const & slimm_db)
{
    std::vector<std::set<uint32_t> > taxa_rank_set;
    taxa_rank_set.resize(LINAGE_LENGTH);

    for(auto ac__taxid_it=slimm_db.ac__taxid.begin(); ac__taxid_it != slimm_db.ac__taxid.end(); ++ac__taxid_it)
    {
        uint32_t tid = ac__taxid_it->second[0];
        if(taxon_ids.find(tid) != taxon_ids.end())
        {
            for (uint32_t i=0; i<LINAGE_LENGTH; ++i)
            {
                taxa_rank_set[i].insert(ac__taxid_it->second[i]);
            }
        }
    }
    for (uint32_t i=0; i<LINAGE_LENGTH; ++i)
    {
        if (taxa_rank_set[i].size() == 1)
            return *(taxa_rank_set[i].begin());
    }
    return 1;
}

uint32_t get_lca(std::set<uint32_t> const & taxon_ids, slimm_database const & slimm_db)
{
    return get_lca(taxon_ids, taxon_ids, slimm_db);
}

uint32_t get_lca(std::vector<uint32_t> const & taxon_ids, slimm_database const & slimm_db)
{
    std::set<uint32_t> s(taxon_ids.begin(), taxon_ids.end());
    if (s.size() == 1)
        return taxon_ids[0];
    else
        return get_lca(s, s, slimm_db);
}


// try to open sam file
inline bool read_bam_file(BamFileIn & bam_file, BamHeader & bam_header, std::string const & bam_file_path)
{
    if (!open(bam_file, toCString(bam_file_path)))
    {
        std::cerr << "Could not open " << bam_file_path << "!\n";
        return false;
    }
    readHeader(bam_header, bam_file);
    return true;
}

inline uint32_t get_avg_read_length(BamFileIn & bam_file, uint32_t const sample_size)
{
    BamAlignmentRecord record;
    uint32_t count = 0, totlaLength = 0;
    while (!atEnd(bam_file) && count < sample_size)
    {
        readRecord(record, bam_file);
        if (length(record.seq) == 0)
            continue;  // Skip records without sequences.
        totlaLength += length(record.seq);
        ++count;
    }
    return totlaLength/count;
}

inline uint32_t get_taxon_id_pos(CharString const & accession)
{
    uint32_t taxon_id_pos = 0;
    if (!get_taxon_id(taxon_id_pos, accession, "ti"))
    {
        if (!get_taxon_id(taxon_id_pos, accession, "kraken:taxid"))
        {
            std::cerr<<"Unable to find a way to resolve taxon id associated with references.\n"
            <<"Make sure you used a set of references provided with SLIMM\n"
            <<"or generated by the preprocessing script.\n";
            exit(1);
        }
    }
    return taxon_id_pos;
}



std::vector<std::string> get_bam_files_in_directory(std::string directory)
{
    std::vector<std::string>  input_paths;
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
           input_paths.push_back(full_file_name);
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
            input_paths.push_back(full_file_name);
        //appendValue(input_paths, full_file_name);
    }
    closedir(dir);
#endif
     return input_paths;
} // get_bam_files_in_directory

