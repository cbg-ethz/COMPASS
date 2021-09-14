#ifndef DEF_NODE
#define DEF_NODE

#include <vector>
#include <set>
#include <string>

#include "Structures.h"
#include "Scores.h"

class Node {
    private: 
        std::vector<int> mutations; //somatic mutations always transform from ref to alt. 

        std::multiset<std::pair<int,std::vector<int>>> CNLOH_events; // pairs (region_index, lost_alleles) 
        // lost_alleles is a vector, where a 0 (resp 1) at position i indicates that at the i-th locus in that region, 
        // a copy of the ref (resp alt) allele is lost and a copy of the other allele is gained
        // I use an (ordered) set, with the default (lexicographic) order, which here sorts the events by position.

        std::multiset<std::tuple<int,int,std::vector<int>>> CNV_events; // triplet (region_index,gain/loss,affected alleles)
        // the second element is +1 in case of copy number gain, -1 in case of copy number loss
        // the third element is a vector of length the number of loci in this region (potentially 0), 
        // and it contains 0 if the ref allele is deleted/duplicated, or 1 for the alt allele
        
        std::set<int> affected_loci; // loci for which the copy number of the ref or alt alleles is changed in the node.
        std::set<int> affected_regions; // regions affected by a CNV event in the node
        bool all_CNV_events_valid; // true if all the CNV events are valid (copy numbers remain >0 etc...)
        
        std::vector<int> n_ref_allele; //array of number of copies of the ref allele, for each variable locus
        std::vector<int> n_alt_allele; //array of number of copies of the alt allele, for each variable locus
        std::vector<int> cn_regions; // copy number of each region

        Scores* cache_scores;


    public:
        std::vector<double> attachment_scores;

        Node(Scores* cache); //constructor
        Node(Node& source); //copy constructor
        Node(Node& source1,Node& source2); // create doublet from 2 nodes
        void init_structures();
        ~Node();

        // Compute genotype of the node, depending on the genotype of its parent
        void update_genotype(Node* parent);
        // compute attachment scores from scratch
        void compute_attachment_scores(bool use_CNV,const std::vector<double>& dropout_rates_ref,
                                        const std::vector<double>& dropout_rates_alt,const std::vector<double>& region_probabilities);
        // just compute difference with parent
        void compute_attachment_scores_parent(bool use_CNV,Node* parent,const std::vector<double>& dropout_rates_ref,
                                        const std::vector<double>& dropout_rates_alt,const std::vector<double>& region_probabilities);

        // MCMC moves for the nodes
        void add_mutation(int locus){mutations.push_back(locus);}
        int remove_random_mutation(); // removes a random mutation, and return the index of the mutation
        void add_CNLOH(std::pair<int,std::vector<int>> CNLOH){CNLOH_events.insert(CNLOH);}
        std::pair<int,std::vector<int>> remove_random_CNLOH();
        void add_CNV(std::tuple<int,int,std::vector<int>> CNV){CNV_events.insert(CNV);}
        std::tuple<int,int,std::vector<int>> remove_random_CNV();
        void remove_CNV(std::tuple<int,int,std::vector<int>> CNV) {CNV_events.erase(CNV_events.lower_bound(CNV));}
        void remove_CNVs_in_region(int region);
        double exchange_CNV_CNLOH(std::vector<int> candidate_regions);
        void change_alleles_CNV();
        
        // Various accessors
        bool is_empty(){ return (mutations.size()==0 && CNLOH_events.size()==0 && CNV_events.size()==0);}

        //  Accessors for Mutations
        int get_n_ref_allele(int locus){return n_ref_allele[locus];}
        int get_n_alt_allele(int locus){return n_alt_allele[locus];}
        std::vector<int> get_mutations() {return mutations;}
        int get_number_mutations() {return mutations.size();}
        
        //  Accessors for CNLOH
        std::multiset<std::pair<int,std::vector<int>>> get_CNLOH_events(){return CNLOH_events;}
        int get_number_CNLOH() {return CNLOH_events.size();}
        int get_number_disjoint_CNLOH(); // here, count contiguous CNLOH as just one event

        //  Accessors for CNVs
        int get_cn_region(int region){return cn_regions[region];}
        std::multiset<std::tuple<int,int,std::vector<int>>> get_CNV_events() {return CNV_events;}
        std::set<int> get_affected_regions(){return affected_regions;}
        bool get_all_CNV_valid(){return all_CNV_events_valid;}
        int get_number_CNV() {return CNV_events.size();}
        int get_number_CNV_mut() {
            int count=0;
            for (auto CNV: CNV_events){
                if (std::get<2>(CNV).size()>0) count+=1;
            }
            return count;
        }
        int get_number_CN_losses(){
            int n_CN_losses = 0;
            for (auto CNV: CNV_events){
                if (std::get<1>(CNV)<0) n_CN_losses++; 
            }
            return n_CN_losses;
        }
        int get_number_disjoint_CNV(std::vector<int> regions_successor);
        int get_number_disjoint_CNV_LOH(std::vector<int> regions_successor); // only count copy number loss of regions containing a variant
        
        
        // Label of the node, used for plotting.
        std::string get_label();
        std::string get_label_simple(std::set<int> excluded_mutations);

};

#endif