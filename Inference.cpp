#include <random>
#include <cfloat>
#include <iostream>
#include <string>

#include "Inference.h"
#include "Structures.h"

#include <assert.h>

extern int n_loci;
extern int n_regions;
extern Data data;
extern Params parameters;

Inference::Inference(std::string name, double temperature):
        cache_scores{new Scores()},
        t{cache_scores,false},
        t_prime{t},
        best_tree{t},
        tree_name{name},
        max_temperature{temperature}
    {

    // Weight of each MCMC move (does not have to be normalized)
    move_weights = {
        1,      // Prune and reattach
        0.5,    // Swap node labels
        1,      // Move mutation
        2,      // Split/merge node
        1,      // Add/remove CNLOH
        1,      // Move CNLOH
        2,      // Add/Remove CNV
        1,      // Move CNV
        0.5,    // Merge or duplicate CNV
        1,      // Exchange CNV/CNLOH
        0.5     // Change alleles affected by CNV
    };
    allow_diff_dropoutrates = true;
}

Inference::~Inference(){
    delete cache_scores;
}

Tree Inference::find_best_tree(bool use_CNV, int nb_steps, int burn_in){

    //First, find the best tree without CNV.
    mcmc(false, nb_steps,burn_in);
    if (!use_CNV){
        if (tree_name!="") best_tree.to_dot(tree_name+".gv");
        return best_tree;
    }

    best_tree.select_regions(); 
    if (!best_tree.contains_candidate_regions()){
        //If cannot find candidate regions which might contain a CNV (or if not cells are attached to the root), return now
        if (tree_name!="") best_tree.to_dot(tree_name+".gv");
        return best_tree;
    }
    

    if (tree_name!="") best_tree.to_dot(tree_name+"_noCNV.gv");
    // Find best tree with CNV
    best_tree.allow_CNV();
    t = best_tree;
    t_prime = t;
    mcmc(true, nb_steps,0);
    if (tree_name!="") best_tree.to_dot(tree_name+".gv");


    return best_tree;
}

void Inference::mcmc(bool use_CNV, int nb_steps,int burn_in){
    best_tree = t;
    double best_log_score = -DBL_MAX;
    int move_id;
    int n_best_tree=0;
    for (int step=0;step<nb_steps;step++){
        if (parameters.verbose) std::cout<<"MCMC step " <<step<<"  ----------------------------------------"<<std::endl;

        int max_move_index=4;
        if (step>=burn_in) max_move_index=6;
        if (step>burn_in && use_CNV) max_move_index = move_weights.size();

        move_id = select_move(max_move_index);
        switch(move_id){
            case 0:
                if (parameters.verbose) std::cout<<"Selected prune and reattach"<<std::endl;
                t_prime.prune_reattach();
                break;
            case 1:
                if (parameters.verbose) std::cout<<"Selected swap node labels"<<std::endl;
                t_prime.swap_node_labels();
                break;
            case 2:
                if (parameters.verbose) std::cout<<"Selected move mutation"<<std::endl;
                t_prime.move_mutation();
                break;
            case 3:
                if (parameters.verbose) std::cout<<"Selected split/merge node"<<std::endl;
                t_prime.split_merge_node();
                break;
            case 4:
                if (parameters.verbose) std::cout<<"Selected add/remove CNLOH"<<std::endl;
                t_prime.add_remove_CNLOH();
                break;
            case 5:
                if (parameters.verbose) std::cout<<"Selected move CNLOH"<<std::endl;
                t_prime.move_CNLOH();
                break;
            case 6:
                if (parameters.verbose) std::cout<<"Selected add/remove CNV"<<std::endl;
                t_prime.add_remove_CNV();
                break;
            case 7:
                if (parameters.verbose) std::cout<<"Selected move CNV"<<std::endl;
                t_prime.move_CNV();
                break;
            case 8:
                if (parameters.verbose) std::cout<<"Selected merge or duplicate CNV"<<std::endl;
                t_prime.merge_or_duplicate_CNV();
                break;
            case 9:
                if (parameters.verbose) std::cout<<"Selected exchange CNV/CNLOH"<<std::endl;
                t_prime.exchange_CNV_CNLOH();
                break;
            case 10:
                if (parameters.verbose) std::cout<<"Selected change alleles affected by CNV"<<std::endl;
                t_prime.change_alleles_CNV();
                break;
        }
        double acceptance_ratio = 0.0;
        if (t_prime.hastings_ratio>0.0){
            t_prime.compute_likelihood(allow_diff_dropoutrates && step>burn_in);
            t_prime.compute_prior_score();
            t_prime.update_full_score();
            double temperature = max_temperature - (max_temperature-1.0) * step / nb_steps;
            acceptance_ratio = std::exp((t_prime.log_score - t.log_score)/temperature + std::log(t_prime.hastings_ratio) /10.0);
        }
        
        double rd = ( (double)std::rand() ) /RAND_MAX;
        if (parameters.verbose)  std::cout<<"t_prime score: " <<t_prime.log_score<<", t score: " << t.log_score
        <<", hastings ratio:"<<t_prime.hastings_ratio<< ", acceptance ratio:" << acceptance_ratio<< ",priorT "<<t.log_prior_score<<", priorT' "<<t_prime.log_prior_score<<std::endl;
        
        if (rd<=acceptance_ratio){
            // accept t_prime
            if (parameters.verbose) std::cout<<"Accepted move" <<std::endl;
            t = t_prime; 
        }
        else{
            // reject t_prime
            if (parameters.verbose) std::cout<<"Rejected move" <<std::endl;
            t_prime = t;
        }
        
        // compare score with best score
        if (t.log_score>best_log_score){
            best_log_score = t.log_score;
            best_tree = t;
            n_best_tree++;
        }
        if (parameters.verbose) std::cout<<"Current best score " <<best_log_score<<std::endl;

    }
}


int Inference::select_move(int max_move_index){
    // Sample a MCMC move, taking into account the weight of each move.
    double sum_weights=0;
    for (int i=0;i<max_move_index;i++) sum_weights+=move_weights[i];
    double rd = sum_weights * ( (double)std::rand() ) /RAND_MAX;
    double cumulative_weight = 0;

    for (int i=0;i<max_move_index;i++){
        cumulative_weight+= move_weights[i];
        if (rd<=cumulative_weight){
            return i;
        }
    }
    return max_move_index;
}