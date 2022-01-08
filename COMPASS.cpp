#include <iostream>
#include <cmath>
#include <cfloat>
#include <random>
#include <fstream>
#include <string.h>
#include <omp.h>

#include "Inference.h"
#include "Tree.h"
#include "Scores.h"
#include "input.h"

int n_cells;
int n_loci;
int n_regions;
std::vector<Cell> cells;
Data data;
Params parameters;

int main(int argc, char* argv[]){
    init_params();
    parameters.verbose=false;
    // Read command line arguments
    std::string input_file{};
    int n_chains=4;
    int chain_length=5000;
    int burn_in = 1000;
    double temperature=10;
    double betabin_overdisp = parameters.omega_het;
    bool use_CNV=true;
    bool apply_filter_regions = true;
    bool output_simplified = true;
    std::string output{};
    data.sex = "female";
    //parameters.verbose=true;
    for (int i=1;i<argc-1;i++){
        if (strcmp(argv[i],"-i")==0){
            input_file = argv[i+1];
        }
        else if (strcmp(argv[i],"--nchains")==0){
            n_chains=atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"--chainlength")==0){
            chain_length=atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"--burnin")==0){
            burn_in=atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"--temperature")==0){
            temperature=atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"--overdisp")==0){
            betabin_overdisp=atof(argv[i+1]);
        }
        else if (strcmp(argv[i],"-o")==0){
            output=argv[i+1];
        }
        else if (strcmp(argv[i],"-d")==0){
            if (strcmp(argv[i+1],"0")==0) parameters.use_doublets=false;
        }
        else if (strcmp(argv[i],"--CNV")==0){
            if (strcmp(argv[i+1],"0")==0) use_CNV=false;
        }
        else if (strcmp(argv[i],"--filterregions")==0){
            if (strcmp(argv[i+1],"0")==0) apply_filter_regions=false;
        }
        else if (strcmp(argv[i],"--sex")==0){
           data.sex= std::string(argv[i+1]);
        }
        else if (strcmp(argv[i],"--prettyplot")==0){
            if (strcmp(argv[i+1],"0")==0) output_simplified=false;
        }
    }

    load_CSV(input_file,use_CNV,apply_filter_regions); 

    parameters.omega_het = std::min(parameters.omega_het,betabin_overdisp);
    parameters.omega_het_indel = std::min(parameters.omega_het_indel,betabin_overdisp);

    // Get the name of the file, without directory
    std::string input_name = input_file;
    int name_start=0;
    for (int i=0;i<input_file.size();i++){
        if (input_file[i]=='/') name_start=i+1;
    }

    input_name=input_file.substr(name_start,input_file.size()-name_start);

	std::vector<double> results{};
    results.resize(n_chains);
    std::vector<Tree> best_trees{};
    best_trees.resize(n_chains);
    if (n_chains<omp_get_num_procs()) omp_set_num_threads(n_chains);
    else omp_set_num_threads(omp_get_num_procs());

    std::cout<<"Starting "<<std::to_string(n_chains)<< " MCMC chains in parallel"<<std::endl;
    #pragma omp parallel for
	for (int i=0;i<n_chains;i++){
		std::srand(i);
		Inference infer{"",temperature,i};
        best_trees[i] = infer.find_best_tree(use_CNV,chain_length,burn_in);
		results[i]=best_trees[i].log_score;
	}
    double best_score=-DBL_MAX;
    int best_score_index=-1;
	for (int i=0;i<n_chains;i++){
		if (best_score<results[i]){
            best_score=results[i];
            best_score_index = i;
        }
	}
    if (output_simplified) best_trees[best_score_index].to_dot_pretty(output);
    else best_trees[best_score_index].to_dot(output);
    std::cout<<"Completed ! The output was written to "<<output<< " in dot format. You can visualize it by running: dot -Tpng "<<output<<" -o output.png"<<std::endl;

	return 0;
}
