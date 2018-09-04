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

#ifndef REFERENCE_CONTIG_H
#define REFERENCE_CONTIG_H


using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================
// ----------------------------------------------------------------------------
// Class taxon_properties
// ----------------------------------------------------------------------------
//class taxon_properties
//{
//    std::string         name;
//    unsigned int        length;
//    unsigned int        level1_neighbors;
//    unsigned int        level2_neighbors;
//    unsigned int        level3_neighbors;
//    float               GCContent;
//
//    taxon_properties(): name(""),
//    length(0),
//    level1_neighbors(0),
//    level2_neighbors(0),
//    level3_neighbors(0),
//    GCContent(50.0) {}
//};

// ----------------------------------------------------------------------------
// Class bins_coverage
// ----------------------------------------------------------------------------
class bins_coverage
{
public:
    uint32_t                    bin_width;
    uint32_t                    number_of_bins;
    std::vector <uint32_t>      bins_height;

    bins_coverage(): bin_width(0),
                     number_of_bins(0) {}

    bins_coverage(uint32_t totalLen, uint32_t width)
    {
        bin_width = width;
        number_of_bins = totalLen/width + 1;
        bins_height.resize(number_of_bins, 0);
    }

    uint32_t none_zero_bin_count()
    {
        if (_none_zero_bin_count == -1)
        {
            _none_zero_bin_count = (number_of_bins - std::count (bins_height.begin(), bins_height.end(), 0));
        }
        return _none_zero_bin_count;
    }

private:
    int32_t _none_zero_bin_count = -1;
};

// ----------------------------------------------------------------------------
// Class reference_contig
// ----------------------------------------------------------------------------
class reference_contig
{
public:
    std::string         accession;
    uint32_t            taxa_id;
    uint32_t            length;
    uint32_t            reads_count;
    uint32_t            uniq_reads_count;
    uint32_t            uniq_reads_count2;
    bins_coverage       cov;
    bins_coverage       uniq_cov;
    bins_coverage       uniq_cov2;
    float               abundance;
    float               uniq_abundance;
    float               uniq_abundance2;

    reference_contig(): accession(""),
                        taxa_id(0),
                        length(0),
                        reads_count(0),
                        uniq_reads_count(0),
                        uniq_reads_count2(0),
                        cov(),
                        uniq_cov(),
                        uniq_cov2(),
                        abundance(0.0),
                        uniq_abundance(0.0),
                        uniq_abundance2(0.0){}

    reference_contig(std::string & ref_name, uint32_t & t_id, uint32_t & ref_length, uint32_t & bin_width):
                        accession(ref_name),
                        taxa_id(t_id),
                        length(ref_length),
                        reads_count(0),
                        uniq_reads_count(0),
                        uniq_reads_count2(0),
                        abundance(0.0),
                        uniq_abundance(0.0),
                        uniq_abundance2(0.0) 
                        {
                            // Intialize coverages based on the length of a refSeq
                            bins_coverage tmp_cov(ref_length, bin_width);
                            cov = tmp_cov;
                            uniq_cov = tmp_cov;
                            uniq_cov2 = tmp_cov;
                        }

    //Member functions
    inline float cov_percent()
    {
        return float(cov.none_zero_bin_count())/cov.number_of_bins;        
    }
    inline float uniq_cov_percent()
    {
        return float(uniq_cov.none_zero_bin_count())/uniq_cov.number_of_bins;        
    }
    inline float uniq_cov_percent2()
    {
        return float(uniq_cov2.none_zero_bin_count())/uniq_cov2.number_of_bins;        
    }

    inline float cov_depth()
    {
        if (_cov_depth == -1)
            _cov_depth = _get_cov_depth(cov);
        return _cov_depth;
    }
    inline float               uniq_cov_depth()
    {
        if (_uniq_cov_depth == -1)
            _uniq_cov_depth = _get_cov_depth(uniq_cov);
        return _uniq_cov_depth;
    }
    inline float uniq_cov_depth2()
    {
        if (_uniq_cov_depth2 == -1)
            _uniq_cov_depth2 = _get_cov_depth(uniq_cov2);
        return _uniq_cov_depth2;
    }

private:
    float _cov_percent = -1;
    float _uniq_cov_percent = -1;
    float _uniq_cov_percent2 = -1;
    float _cov_depth = -1;
    float _uniq_cov_depth = -1;
    float _uniq_cov_depth2 = -1;

    // --------------------------------------------------------------------------
    // Function getCovDepth()
    // --------------------------------------------------------------------------
    inline float _get_cov_depth(bins_coverage & c)
    {
        if(c.none_zero_bin_count() == 0)
            return 0.0;

        std::vector <uint32_t>::iterator it;

        //copy the coverage height.
        std::vector <float>  bHeights;

        bHeights.reserve(c.number_of_bins);
        for (unsigned int i=0; i<c.number_of_bins; ++i )
        {
            bHeights.push_back(float(c.bins_height[i]));
        }
        return mean(bHeights);
    }
};
#endif /* REFERENCE_CONTIG_H */
