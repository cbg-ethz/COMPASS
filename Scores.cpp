#include <cmath>
#include <vector>
#include <cfloat>
#include <string>
#include <iostream>
#include <random>

#include "Scores.h"
#include "Node.h"
#include "Structures.h"


//global variables
extern int n_loci;
extern int n_cells;
extern int n_regions;

extern std::vector<Cell> cells;
extern Data data;
extern Params parameters;



Scores::Scores(){
    std::map<long int,std::vector<std::vector<double>>> cache_likelihood_allelecounts_cells{};
    std::map<int,std::map<long int,std::vector<double>>> cache_dropoutrate_cellscores{};
    std::map<long int,std::vector<double>> cache_dropoutscores{};

    std::map<int,std::map<long int,std::vector<double>>> cache_dropoutrate_dropoutsref{};
    std::map<int,std::map<long int,std::vector<double>>> cache_dropoutrate_dropoutsalt{};
    std::map<int,double> cache_n_choose_k{};
    count_cache=0;
}

double Scores::log_sum_exp(const std::vector<double>& terms){
    double max = -DBL_MAX;
    for (int i=0;i<terms.size();i++){
        if (terms[i]>max) max=terms[i];
    }

    double sum=0;
    for (int i=0;i<terms.size();i++){
        sum+= std::exp(terms[i]-max);
    }
    return std::log(sum) + max;
}

std::vector<double> Scores::log_sum_exp_vector(const std::vector<std::vector<double>>& terms){
    std::vector<double> maxs(n_cells,-DBL_MAX);
    for (int i=0;i<terms.size();i++){
    	for (int j=0;j<n_cells;j++){
    		if (terms[i][j]>maxs[j]) maxs[j]=terms[i][j];
    	}
    }

    std::vector<double> sums(n_cells,0.0);
    for (int i=0;i<terms.size();i++){
    	for (int j=0;j<n_cells;j++){
    		sums[j]+= std::exp(terms[i][j]-maxs[j]);
    	}
    }
    for (int j=0;j<n_cells;j++) sums[j] = std::log(sums[j]) + maxs[j];
    return sums;
}

double Scores::log_n_choose_k(int n, int k){
    int hash = n+1000*k;
    if (cache_n_choose_k.count(hash)) return cache_n_choose_k[hash];
    double res=0;
    for (int i=0; i<k;i++){
        res+= std::log(n-i);
        res-= std::log(i+1);
    }
    cache_n_choose_k[hash] = res;
    return res;
}


std::vector<double> Scores::compute_SNV_loglikelihoods(int c_ref,int c_alt,int locus, double dropout_rate_ref, double dropout_rate_alt){
    // compute the SNV log-likelihood for one locus for all of the cells
    // c_ref and c_alt are the copy numbers of each allele in the genotype
    
    // If homozygous, the copy number of the only allele is irrelevant for the allelic proportion
    if (c_ref==0) c_alt=1;
    else if (c_alt==0) c_ref==1;

    long int hash = c_ref+ 20*c_alt + 400*locus ;
    
    int discretized_dropout_rate_ref = std::round(dropout_rate_ref*1000);
    int discretized_dropout_rate_alt = std::round(dropout_rate_alt*1000);
    int hash_dropoutrates = discretized_dropout_rate_ref + 1000 * discretized_dropout_rate_alt;
    // if we already computed the scores for this dropout rate and this combination of locus, c_ref and c_alt.
    if (cache_dropoutrate_cellscores.count(hash_dropoutrates) && cache_dropoutrate_cellscores[hash_dropoutrates].count(hash)){
        return cache_dropoutrate_cellscores[hash_dropoutrates][hash];
    }
    dropout_rate_ref = 1.0*discretized_dropout_rate_ref/1000.0;
    dropout_rate_alt = 1.0*discretized_dropout_rate_alt/1000.0;
    std::vector<std::vector<double>> likelihood_alleles_cells{};
    likelihood_alleles_cells.resize((c_ref+1)*(c_alt+1)-1);
    std::vector<double> dropoutscores{};
    dropoutscores.resize(n_cells);
    int idx=0;
    for (int k=0;k<=c_ref;k++){
        for (int l=0;l<=c_alt;l++){
            if (k==0 & l==0) continue;
            long int hash_dropout = k+20*l + 400 * locus;
            if (cache_dropoutscores.count(hash_dropout)){
                dropoutscores = cache_dropoutscores[hash_dropout];                
            }
            else{ // does not depend on dropout rate so does not have to be computed often
                std::vector<double> seq_error_rates(n_cells,0.0);
                double f;
                double likelihood_dropout;
                double omega;
                double eps1; // sequencing errors from ref to alt (always rare)
                double eps2; // sequencing errors from alt to ref (seem quite common for indels)
                if (data.variant_is_SNV[locus]){
                    eps1 = parameters.sequencing_error_rate;
                    eps2 = parameters.sequencing_error_rate;
                    if (k==0 || l==0) omega = parameters.omega_hom;
                    else omega = parameters.omega_het;
                }
                else{
                    eps1 = parameters.sequencing_error_rate;
                    eps2 = parameters.sequencing_error_rate_indel;
                    if (k==0 || l==0) omega = parameters.omega_hom_indel;
                    else omega = parameters.omega_het_indel;
                }

                f = 1.0*l/(k+l) * (1-eps2) + 1.0*k/(k+l) * eps1; // frequency of the alt nucleotide
                for (int j=0;j<n_cells;j++){
                    int ref_count= cells[j].ref_counts[locus];
                    int alt_count=cells[j].alt_counts[locus];
                    likelihood_dropout = std::lgamma(alt_count + omega*f) + std::lgamma(ref_count + omega*(1-f)) - std::lgamma(alt_count+ref_count + omega)
                    -std::lgamma(omega*f) - std::lgamma(omega*(1-f)) + std::lgamma(omega);
                    //the term (a+r) choose a does not depend on the parameters, so does not need to be included here. + log_n_choose_k(ref_count+alt_count,alt_count)

                    dropoutscores[j] = likelihood_dropout;
                }
                cache_dropoutscores[hash_dropout] = dropoutscores;
            }
            // add dropout probability
            likelihood_alleles_cells[idx].resize(n_cells);
            double dropout_prob = log_n_choose_k(c_ref,k) + log_n_choose_k(c_alt,l) 
                                            +(c_ref-k) * std::log(dropout_rate_ref) + k * std::log(1-dropout_rate_ref)
                                            + (c_alt -l)*std::log(dropout_rate_alt) + l * std::log(1-dropout_rate_alt);
            // Cannot have a dropout of all the alleles
            double all_dropout_prob = std::pow(dropout_rate_ref,c_ref) * std::pow(dropout_rate_alt,c_alt);
            dropout_prob-= std::log(1-all_dropout_prob);
            for (int j=0;j<n_cells;j++){
                likelihood_alleles_cells[idx][j] = dropoutscores[j] + dropout_prob;
            }
            idx+=1;
        }
    }
    // Compute the score for each cell by summing over each dropout combination
    std::vector<double> scores_cells = log_sum_exp_vector(likelihood_alleles_cells);
    std::vector<double> dropoutsref(n_cells,0.0);
    std::vector<double> dropoutsalt(n_cells,0.0);
    double config_prob;

    // Compute the of average number of dropouts at this locus that occured in each cell
    // This is used for optimizing the dropout rate
    idx=0;
    for (int k=0;k<=c_ref;k++){
        for (int l=0;l<=c_alt;l++){
            if (k==0 && l==0) continue;
            for (int j=0;j<n_cells;j++){
                config_prob = std::exp(likelihood_alleles_cells[idx][j]-scores_cells[j]);
                dropoutsref[j]+=1.0*(c_ref-k) * config_prob;
                dropoutsalt[j]+=1.0*(c_alt-l) * config_prob;
            }
            idx+=1;
        }
    }

    cache_dropoutrate_cellscores[hash_dropoutrates][hash] = scores_cells;
    cache_dropoutrate_dropoutsref[hash_dropoutrates][hash] = dropoutsref;
    cache_dropoutrate_dropoutsalt[hash_dropoutrates][hash] = dropoutsalt;
    count_cache++;
    return scores_cells;
}

std::vector<double> Scores::compute_CNV_loglikelihoods(int region, double region_proportion){
    // Compute the likelihood of the read count in the region, based on the negative binomial distribution (Gamma-Poisson)
    // theta is the scale parameter for the Gamma distribution.
    int discretized_region_proportion = std::round(region_proportion*1000);
    long int hash = region + 500*discretized_region_proportion;
    if (cache_cnvlikelihood_cells.count(hash)) return cache_cnvlikelihood_cells[hash];


    double theta = parameters.theta;
    std::vector<double> cnv_loglikelihoods{};
    cnv_loglikelihoods.resize(n_cells);
    for (int j=0;j<n_cells;j++){
        double expected_read_count_region = cells[j].total_counts * region_proportion;
        cnv_loglikelihoods[j] = std::lgamma(cells[j].region_counts[region] + theta-1) + theta * std::log(theta / (theta + expected_read_count_region))
                                + cells[j].region_counts[region] * std::log(expected_read_count_region / (expected_read_count_region+theta));
    }
    cache_cnvlikelihood_cells[hash] = cnv_loglikelihoods;
    count_cache++;
    return cnv_loglikelihoods;
}



std::vector<double> Scores::get_dropoutref_counts_genotype(int c_ref,int c_alt, int locus,double dropout_rate_ref,double dropout_rate_alt){
    if (c_ref==0) c_alt=1;
    else if (c_alt==0) c_ref==1;
    if (c_ref==0 || c_alt==0){
        dropout_rate_ref=0.1;
        dropout_rate_alt=0.1;
    }
    long int hash = c_ref+ 20*c_alt + 400*locus;

    int discretized_dropout_rate_ref = std::round(dropout_rate_ref*1000);
    int discretized_dropout_rate_alt = std::round(dropout_rate_alt*1000);
    int hash_dropoutrate = discretized_dropout_rate_ref + discretized_dropout_rate_alt * 1000;
    // if we already computed the scores for this dropout rate and this combination of locus, c_ref and c_alt.
    if (cache_dropoutrate_dropoutsref.count(hash_dropoutrate) && cache_dropoutrate_dropoutsref[hash_dropoutrate].count(hash)){
        return cache_dropoutrate_dropoutsref[hash_dropoutrate][hash];
    }
    else{
        auto temp = compute_SNV_loglikelihoods(c_ref,c_alt,locus,dropout_rate_ref,dropout_rate_alt);
        return cache_dropoutrate_dropoutsref[hash_dropoutrate][hash];
    }
}

std::vector<double> Scores::get_dropoutalt_counts_genotype(int c_ref,int c_alt, int locus,double dropout_rate_ref,double dropout_rate_alt){
    if (c_ref==0) c_alt=1;
    else if (c_alt==0) c_ref==1;
    if (c_ref==0 || c_alt==0){
        dropout_rate_ref=0.1;
        dropout_rate_alt=0.1;
    }
    long int hash = c_ref+ 20*c_alt + 400*locus ;

    int discretized_dropout_rate_ref = std::round(dropout_rate_ref*1000);
    int discretized_dropout_rate_alt = std::round(dropout_rate_alt*1000);
    int hash_dropoutrate = discretized_dropout_rate_ref + discretized_dropout_rate_alt * 1000;
    // if we already computed the scores for this dropout rate and this combination of locus, c_ref and c_alt.
    if (cache_dropoutrate_dropoutsalt.count(hash_dropoutrate) && cache_dropoutrate_dropoutsalt[hash_dropoutrate].count(hash)){
        return cache_dropoutrate_dropoutsalt[hash_dropoutrate][hash];
    }
    else{
        auto temp = compute_SNV_loglikelihoods(c_ref,c_alt,locus,dropout_rate_ref,dropout_rate_alt);
        return cache_dropoutrate_dropoutsalt[hash_dropoutrate][hash];
    }
}


void Scores::clear_cache(){
    cache_dropoutrate_cellscores.clear();
    cache_dropoutrate_dropoutsalt.clear();
    cache_dropoutrate_dropoutsref.clear();

    cache_cnvlikelihood_cells.clear();

    count_cache=0;
}

void Scores::clear_cache_if_too_large(){
    // To avoid using too much memory
    if (count_cache>8000) clear_cache();
}





