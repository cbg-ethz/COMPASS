#ifndef DEF_INFERENCE
#define DEF_INFERENCE

#include <vector>
#include <string>
#include "Tree.h"

class Inference{
    private:
        Scores* cache_scores;
        Tree t; // current tree
        Tree t_prime; // proposed tree
        Tree best_tree; // best tree seen so far
        bool allow_diff_dropoutrates;
        double max_temperature; // start MCMC at this temperature, and decreases to 1
        std::vector<double> move_weights;
        std::string tree_name;
        
    public:
        Inference(std::string name="best_tree.gv", double temperature=10.0);
        ~Inference();
        Tree find_best_tree(bool use_CNV=true,int nb_steps=5000, int burn_in=1000);
        void mcmc(bool use_CNV, int nb_steps, int burn_in=1000);
        int select_move(int max_move_index);

};

#endif