#ifndef DEF_SCORES
#define DEF_SCORES

#include <vector>
#include <map>

#include "Structures.h"



class Scores{
    private:
        // Cache various results to avoid computing them repeatedly
        std::map<int,double> cache_n_choose_k;
        std::map<long int,std::vector<std::vector<double>>> cache_likelihood_allelecounts_cells;
        std::map<long int,std::vector<double>> cache_dropoutscores;
        std::map<int,std::map<long int,std::vector<double>>> cache_dropoutrate_cellscores;
        std::map<int,std::map<long int,std::vector<double>>> cache_dropoutrate_dropoutsref;
        std::map<int,std::map<long int,std::vector<double>>> cache_dropoutrate_dropoutsalt;
        std::map<long int, std::vector<double>> cache_cnalikelihood_cells;
        int count_cache;
    public:
        Scores();
        double log_n_choose_k(int n, int k);
        static double log_sum_exp(const std::vector<double>& terms); //performs a robust log_sum_exp (ie sum exp(x-max) for better precision)
        static std::vector<double> log_sum_exp_vector(const std::vector<std::vector<double>>& terms);
        //static double compute_SNV_loglikelihood(int c_ref,int c_alt, int ref_count, int alt_count, double eps, double omega,double mu);
        std::vector<double> compute_SNV_loglikelihoods(int c_ref,int c_alt, int locus, double dropout_rate_ref, double dropout_rate_alt);


        std::vector<double> compute_CNA_loglikelihoods(int region, double region_probability);

        std::vector<double> get_dropoutref_counts_genotype(int c_ref,int c_alt, int locus,double dropout_rate_ref,double dropout_rate_alt);
        std::vector<double> get_dropoutalt_counts_genotype(int c_ref,int c_alt, int locus,double dropout_rate_ref,double dropout_rate_alt);
        //std::vector<double> get_regiondropout_probs(int region, double region_probability);
        //std::vector<double> get_default_CNV_scores();
        void clear_cache();   
        void clear_cache_if_too_large();
        int get_cache_size(){return count_cache;}

};

#endif