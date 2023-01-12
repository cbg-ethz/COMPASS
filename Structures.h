#ifndef DEF_STRUCTURES
#define DEF_STRUCTURES
#include <vector>
#include <map>


struct Cell{
    std::vector<int> ref_counts; // number of ref reads for each variable locus
    std::vector<int> alt_counts; // number of alt reads for each variable locus
    std::vector<int> region_counts; // number of reads in each region
    std::vector<int> genotypes;
    std::vector<int> GQ;
    std::string name;
    int total_counts; // sum of the read counts in each region
};

struct Data{
    std::string sex; // "male" (1 chr X) or "female" (2 chr X)
    std::vector<std::string> locus_to_chromosome;
    std::vector<int> locus_to_position;
    std::vector<std::string> locus_to_reference;
    std::vector<std::string> locus_to_alternative;
    std::vector<bool> variant_is_SNV; // true if the variant is a SNV, false otherwise (indel)
    std::vector<double> locus_to_freq; // frequency of the variant in the population
    std::vector<std::string> locus_to_name; //name of the variant, which describes the effect on the protein (if any)
    std::vector<int> locus_to_region;
    
    std::vector<std::vector<int>> region_to_loci; //list of loci on this amplicon (possibly empty)
    std::vector<std::string> region_to_chromosome;
    std::vector<std::string> region_to_name;
    std::vector<bool> region_is_reliable; // true for the amplicons that are used for computing the CNV score.

};

struct Params{
    double sequencing_error_rate;
    double omega_hom; //concentration parameter for beta-binomial in the SNV likelihood, when homozygous
    double omega_het; //concentration parameter for beta-binomial in the SNV likelihood, when heterozygous

    double sequencing_error_rate_indel; // more errors for indels than SNVs
    double omega_hom_indel;
    double omega_het_indel;

    double prior_dropoutrate_mean; // beta prior for the dropout rates. the mean is alpha/(alpha+beta)
    double prior_dropoutrate_omega; // beta prior for the dropout rates. the inverse overdispersion parameter omega is alpha+beta

    double theta; // inverse dispersion of the negative binomial for the CNV likelihood

    double doublet_rate;
    
    bool use_doublets;
    bool filter_regions;
    bool filter_regions_CNLOH;
    bool verbose;

    // Penalties in the tree prior
    double node_cost;
    double CNA_cost;
    double LOH_cost;
    double mut_notAtRoot_cost;
    double mut_notAtRoot_freq_cost;
};

#endif