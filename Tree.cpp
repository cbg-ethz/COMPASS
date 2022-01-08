#include <random>
#include <cmath>
#include <cfloat>
#include <vector>
#include <stack>
#include <iostream>
#include <fstream>

#include "Tree.h"
#include "Node.h"
#include "Scores.h"
#include "input.h"

#include <assert.h>


//global variables
extern int n_cells;
extern int n_loci;
extern int n_regions;
extern std::vector<Cell> cells;
extern Data data;
extern Params parameters;

Tree::Tree(Scores* cache, bool use_CNV): 
    hastings_ratio(-1.0),
    use_CNV(use_CNV)
{
    cache_scores = cache;
    n_nodes = 6;
    dropout_rates = std::vector<double>(n_loci,0.05);
    dropout_rates_ref = std::vector<double>(n_loci,0.05);
    dropout_rates_alt = std::vector<double>(n_loci,0.05);

    //Generate a random tree
    parents.resize(n_nodes);
    parents[0]=-1; //root
    for (int i=1;i<n_nodes;i++){
        parents[i]=std::rand()%i; // new node attached randomly to one of the existing nodes.
    }
    children.resize(n_nodes);
    compute_children();

    //Initialize the nodes.
    nodes.resize(n_nodes);
    node_probabilities.resize(n_nodes);
    for (int i=0;i<n_nodes;i++){
        nodes[i] = new Node(cache_scores);
        node_probabilities[i] = 1.0/n_nodes;
    }

    //randomly assign each somatic mutation to a node and compute the genotypes
    for (int i=0;i<n_loci;i++){
        nodes[std::rand()%n_nodes]->add_mutation(i);
    }
    compute_nodes_genotypes();

    // Initialize utils
    cells_attach_loglik.resize(n_cells);
    cells_loglik.resize(n_cells);
    cells_attach_prob.resize(n_cells);
    best_attachments.resize(n_cells);

    compute_likelihood();
    compute_prior_score();
    update_full_score();
}

Tree::Tree(){
    //nodes.clear();
    //doublets.clear();
    cells_attach_loglik.resize(n_cells);
    cells_loglik.resize(n_cells);
    cells_attach_prob.resize(n_cells);
    best_attachments.resize(n_cells);
}


Tree::Tree(const Tree& source):  
    // copy constructor
    use_CNV(source.use_CNV),
    n_nodes(source.n_nodes),
    parents(source.parents),
    children(source.children),
    DFT_order(source.DFT_order),
    node_probabilities(source.node_probabilities),
    dropout_rates(source.dropout_rates),
    dropout_rates_ref(source.dropout_rates_ref),
    dropout_rates_alt(source.dropout_rates_alt),
    region_probabilities(source.region_probabilities),
    candidate_regions(source.candidate_regions),
    regions_successor(source.regions_successor),
    log_prior_score(source.log_prior_score),
    log_likelihood(source.log_likelihood),
    log_score(source.log_score),
    cache_scores(source.cache_scores)
{
    // Deep copy of the nodes
    for (Node* node: source.nodes){
        nodes.push_back(new Node(*node));
    }

    // Initialize utils
    cells_attach_loglik.resize(n_cells);
    cells_loglik.resize(n_cells);
    cells_attach_prob.resize(n_cells);
    best_attachments = source.best_attachments;
}



Tree::~Tree(){
    for (Node* node: nodes){
        delete node;
    }
    for (Node* doublet: doublets){
        delete doublet;
    }
}

Tree& Tree::operator=(const Tree& source){
    // delete the old nodes
    for (Node* n: nodes){
        delete n;
    }
    nodes.clear();
    if (parameters.use_doublets){
        for (Node* doublet: doublets){
            delete doublet;
        }
        doublets.clear();
    }
    use_CNV = source.use_CNV;
    cache_scores = source.cache_scores;
    n_nodes = source.n_nodes;
    parents = source.parents;
    children = source.children;
    DFT_order=source.DFT_order;
    node_probabilities = source.node_probabilities;
    
    dropout_rates = source.dropout_rates;
    dropout_rates_ref = source.dropout_rates_ref;
    dropout_rates_alt = source.dropout_rates_alt;
   
    region_probabilities = source.region_probabilities;
    candidate_regions = source.candidate_regions;
    regions_successor = source.regions_successor;
    
    log_prior_score = source.log_prior_score;
    log_likelihood = source.log_likelihood;
    log_score = source.log_score;
    best_attachments = source.best_attachments;
    
    for (int i=0;i<n_nodes;i++){
        nodes.push_back(new Node(*source.nodes[i]));
    }
    return *this;
}

void Tree::compute_children(){
    // Compute the list of children of each node, from the parent vector.
    children.resize(n_nodes);
    for (int i=0; i <n_nodes; i++){
        children[i].clear();
    }
    for (int i=1; i <n_nodes; i++){ //start from 1 because the root has no parent  
        children[parents[i]].push_back(i);
    }

    // Also compute the DFT order
    std::stack<int> stk;
    DFT_order.clear();
    stk.push(0);
    while (!stk.empty()) {
        int top = stk.top();
        DFT_order.push_back(top);
        stk.pop();
        for (int child: children[top]) {
            stk.push(child);
        };
    }
}

bool Tree::is_ancestor(int potential_ancestor, int potential_descendant){
    if (potential_ancestor==0) return true;
    int ancestor = potential_descendant;
    while (ancestor!=0){
        if (ancestor==potential_ancestor) return true;
        ancestor = parents[ancestor];
    }
    return false;
}


void Tree::compute_nodes_genotypes(){
    // Perform a depth-first traversal and compute the genotype of each node, based on the genotype of its parent.
    std::stack<int> stk;
    stk.push(0); 
    while (!stk.empty()) {
        int top = stk.top();
        stk.pop();
        for (int child: children[top]) {
            stk.push(child);
        };
        if (top!=0) nodes[top]->update_genotype(nodes[parents[top]]);
        else nodes[top]->update_genotype(nullptr);
    }
    if (parameters.use_doublets){
        for (Node* doublet: doublets){
            delete doublet;
        }
        doublets.clear();
        for (int k=0;k<n_nodes;k++){ //k: first node in the doublet
            for (int l=k;l<n_nodes;l++){ //l: second node in the doublet
                // the first node in the doublet should be the one which is earlier is the DFT order
                // when we compute the scores of the doublet from its parent, we only update the updated loci of the second node
                int i=0;
                while (DFT_order[i]!=k && DFT_order[i]!=l) i++;
                if (DFT_order[i]==k) doublets.push_back(new Node(*nodes[k],*nodes[l]));
                else doublets.push_back(new Node(*nodes[l],*nodes[k]));
                
            }
        }
    }
}




void Tree::compute_attachment_scores(bool use_doublets_local){
    // Compute the attachment scores of the nodes below node_index (including node_index) by performing a DFT
    std::stack<int> stk;
    stk.push(0);
    while (!stk.empty()) {
        int top = stk.top();
        stk.pop();
        for (int child: children[top]) {
            stk.push(child);
        };
        if (top==0) // root: compute score from scratch
            nodes[top]->compute_attachment_scores(use_CNV,dropout_rates_ref,dropout_rates_alt,region_probabilities);
        else // start from the parent score, and only compute the difference on loci/regions affected by events
            nodes[top]->compute_attachment_scores_parent(use_CNV,nodes[parents[top]], dropout_rates_ref,dropout_rates_alt,region_probabilities);
    }

    if (use_doublets_local){
        for (int k=0;k<n_nodes;k++){
            for (int l=k;l<n_nodes;l++){ 
                // traverse in DFT order, so that the scores of the parent are always already computed
                int n1 = DFT_order[k]; // first node of the doublet
                int n2 = DFT_order[l]; // second node of the doublet
                int idx; // index in the doublets vector. The vector only contains pairs (k,l) with k<=l
                if (n2>n1) idx = n_nodes * n1 - (n1*(n1-1))/2 + n2-n1;
                else idx = n_nodes * n2 - (n2*(n2-1))/2 + n1-n2;
                // compute doublet (0,0) from scratch
                if (n1==0 && n2==0) doublets[idx]->compute_attachment_scores(use_CNV,dropout_rates_ref,dropout_rates_alt,region_probabilities);
                // compute doublet (n1,n1) from (parent(n1),n1) (or (n1,parent(n1) )
                else if (n1==n2){
                    int idx_parent;
                    if (parents[n1]>n1) idx_parent = n_nodes * n1 - (n1*(n1-1))/2 + parents[n1]-n1;
                    else idx_parent = n_nodes * parents[n1] - (parents[n1]*(parents[n1]-1))/2 + n1- parents[n1];
                    doublets[idx]->compute_attachment_scores_parent(use_CNV,doublets[idx_parent],dropout_rates_ref,dropout_rates_alt,region_probabilities);
                } 
                // compute doublet (n1,n2) from (n1,parent(n2)) (or (parent(n2),n1) )
                else {
                    int idx_parent;
                    if (parents[n2]>n1) idx_parent = n_nodes * n1 - (n1*(n1-1))/2 + parents[n2]-n1;
                    else idx_parent = n_nodes * parents[n2] - (parents[n2]*(parents[n2]-1))/2 + n1- parents[n2];
                    doublets[idx]->compute_attachment_scores_parent(use_CNV,doublets[idx_parent],dropout_rates_ref,dropout_rates_alt,region_probabilities);
                }
            }
        }
    }

}



void Tree::compute_likelihood(bool allow_diff_dropoutrates){
    // Compute the likelihood by marginalizing over the attachment points
    compute_nodes_genotypes();
    cache_scores->clear_cache_if_too_large();

    // EM algorithm to find best node probabilities and dropout rates
    // start with uniform node probabilities and the previous dropout rates
    for (int k=0;k<n_nodes;k++) node_probabilities[k]=1.0/n_nodes;

    avg_diff_nodeprob=10.0; // how much the node probabilities changed between 2 EM steps
    avg_diff_dropoutrates=10.0; // how much the dropout rates changed between 2 EM steps

    bool use_doublets_EM = false;
    int n_loops=0;
    while((avg_diff_nodeprob>0.0005|| avg_diff_dropoutrates>0.0001) && n_loops<100){
        if (avg_diff_dropoutrates>0.0005) compute_attachment_scores(use_doublets_EM); // attachment scores of cells to nodes do not depend on node probabilities
        compute_cells_likelihoods(use_doublets_EM);
        EM_step(use_doublets_EM,false);
        n_loops++;
    }
    // See if likelihood can be improved by allowing, for some loci, the 2 alleles to have different dropout rates
    if (allow_diff_dropoutrates){
        EM_step(use_doublets_EM,true);
        n_loops=0;
        while((avg_diff_nodeprob>0.0005|| avg_diff_dropoutrates>0.0001) && n_loops<40){
            if (avg_diff_dropoutrates>0.0005) compute_attachment_scores(use_doublets_EM); // attachment scores of cells to nodes do not depend on node probabilities
            compute_cells_likelihoods(use_doublets_EM);
            EM_step(use_doublets_EM,true);
            n_loops++;
        }
    }
    if (parameters.use_doublets && !use_doublets_EM){
        // if we did not use doublets for the EM algorithm, we need to recompute the scores with the doublets
        compute_attachment_scores(parameters.use_doublets); 
        compute_cells_likelihoods(parameters.use_doublets);
    }
    log_likelihood=0;
    for (int j=0;j<n_cells;j++){
        log_likelihood+=cells_loglik[j];
    }
}


void Tree::compute_cells_likelihoods(bool use_doublets_local){
    int n_attachment_points = n_nodes;
    if (use_doublets_local) n_attachment_points = (n_nodes*(n_nodes+3)) / 2; //can attach to a node or to a doublet

    for (int j=0; j<n_cells; j++){ 
        double max_loglik=-DBL_MAX;
        int best_attach=-1;
        cells_attach_loglik[j].resize(n_attachment_points);
        for (int k=0; k<n_nodes;k++){
            cells_attach_loglik[j][k] = nodes[k]->attachment_scores[j] + std::log(node_probabilities[k]);
            if (use_doublets_local) cells_attach_loglik[j][k]+=std::log(1-parameters.doublet_rate);
            if (cells_attach_loglik[j][k] > max_loglik){
                max_loglik=cells_attach_loglik[j][k];
                best_attach=k;
            }
        }
        if (use_doublets_local){
            int idx=0;
            for (int k=0;k<n_nodes;k++){
                for (int l=k;l<n_nodes;l++){
                    cells_attach_loglik[j][n_nodes+idx] = doublets[idx]->attachment_scores[j] + std::log(node_probabilities[k]) 
                    + std::log(node_probabilities[l]) +std::log(parameters.doublet_rate);
                    if (k!=l)  cells_attach_loglik[j][n_nodes+idx] += std::log(2); // the doublet (l,k) has the same probability as (k,l)
                    if ( cells_attach_loglik[j][n_nodes+idx] > max_loglik){
                        max_loglik =  cells_attach_loglik[j][n_nodes+idx];
                        best_attach = n_nodes+idx;
                    }
                    idx++;
                }
            }
        }
        cells_loglik[j] = Scores::log_sum_exp(cells_attach_loglik[j]);
        best_attachments[j] = best_attach;
    }
}

void Tree::EM_step(bool use_doublets_local, bool allow_diff_dropoutrates){
    int n_attachment_points = n_nodes;
    if (use_doublets_local) n_attachment_points = (n_nodes*(n_nodes+3)) / 2; //can attach to a node or to a doublet
    nodes_attachment_counts.resize(n_nodes);
    for (int k=0;k<n_nodes;k++) nodes_attachment_counts[k]=0.0;
    //E-step: compute probabilities that cell j is attached to node k, number of cells attached to each node and number of dropouts
    for (int j=0; j<n_cells; j++){ 
        cells_attach_prob[j].resize(n_attachment_points);

        for (int k=0;k<n_nodes;k++){
            cells_attach_prob[j][k] = std::exp(cells_attach_loglik[j][k]-cells_loglik[j]);
            nodes_attachment_counts[k]+=cells_attach_prob[j][k]; // add posterior probability that cell j is attached to node k
        }

        if (use_doublets_local){
            int idx=0;
            for (int k=0;k<n_nodes;k++){
                for (int l=k;l<n_nodes;l++){
                    cells_attach_prob[j][n_nodes+idx] = std::exp(cells_attach_loglik[j][n_nodes+idx]-cells_loglik[j]);
                    nodes_attachment_counts[k] += cells_attach_prob[j][n_nodes+idx]/2.0;
                    nodes_attachment_counts[l] += cells_attach_prob[j][n_nodes+idx]/2.0;
                    idx++;
                }
            }
        }
    }
    // count dropouts 
    double doublet_rate_local = parameters.doublet_rate;
    if (!use_doublets_local) doublet_rate_local=0.0;
    std::vector<double> dropoutref_counts(n_loci,0.0);
    std::vector<double> dropoutalt_counts(n_loci,0.0);
    std::vector<double> nrefalleles(n_loci,0.0);
    std::vector<double> naltalleles(n_loci,0.0);

    for (int i=0;i<n_loci;i++){
        std::map<std::pair<int,int>,std::vector<double>> genotypes_prob; // for each cell, probability that it is attached to a node with a given genotype for locus i
        // First, group nodes depending on their genotype at locus i
        for (int k=0;k<n_nodes;k++){
            int n_ref=nodes[k]->get_n_ref_allele(i);
            int n_alt=nodes[k]->get_n_alt_allele(i);
            std::pair<int,int> p=  std::make_pair(n_ref,n_alt); //p: genotype at locus i for node k
            if (!genotypes_prob.count(p)){
                genotypes_prob[p]=std::vector<double>(n_cells,0.0);
            }
            for (int j=0;j<n_cells;j++){
                genotypes_prob[p][j]+=cells_attach_prob[j][k] / (1.0-doublet_rate_local); // do not use doublets for inferring dropout rates (faster)
            }
        }

        // Count dropouts for each genotype at locus i
        for (const auto& m: genotypes_prob){
            std::pair<int,int> p=m.first;
            int n_ref=p.first;
            int n_alt=p.second;
            if (n_ref>0 && n_alt>0){ // only consider dropouts for heterozygous genotypes. 
                const std::vector<double>& dropoutsref= cache_scores->get_dropoutref_counts_genotype(n_ref,n_alt,i,dropout_rates_ref[i],
                                                                                                    dropout_rates_alt[i]);
                const std::vector<double>& dropoutsalt= cache_scores->get_dropoutalt_counts_genotype(n_ref,n_alt,i,dropout_rates_ref[i],
                                                                                                    dropout_rates_alt[i]);
                for (int j=0;j<n_cells;j++){
                    // if we have no reads at this position in this cell, we can't tell if a dropout occurred.
                    if (cells[j].ref_counts[i]+cells[j].alt_counts[i]>4){
                        dropoutref_counts[i]+=genotypes_prob[p][j] * dropoutsref[j];
                        dropoutalt_counts[i]+=genotypes_prob[p][j] * dropoutsalt[j];
                        nrefalleles[i]+=genotypes_prob[p][j] * n_ref; 
                        naltalleles[i]+=genotypes_prob[p][j] * n_alt; 
                    }
                }
            }
        }   
    }
    

    //M-step: optimize the node probabilities and dropout rates
    avg_diff_nodeprob=0.0;
    for (int k=0;k<n_nodes;k++){ // update node probabilities
        double prev = node_probabilities[k];
        node_probabilities[k] = nodes_attachment_counts[k] / n_cells;
        avg_diff_nodeprob+=std::abs(node_probabilities[k]-prev) / n_nodes;
    } 

    avg_diff_dropoutrates=0.0;
    for (int i=0;i<n_loci;i++){ // update dropout rates
        double prev_ref = dropout_rates_ref[i];
        double prev_alt = dropout_rates_alt[i];
        dropout_rates[i] = ((parameters.prior_dropoutrate_mean*parameters.prior_dropoutrate_omega-1)*2 + dropoutref_counts[i]+dropoutalt_counts[i]) 
                                    / ((parameters.prior_dropoutrate_omega-2)*2+nrefalleles[i]+naltalleles[i]);
        if (allow_diff_dropoutrates){
            // See if likelihood can be improved by allowing 2 different dropout rates
            dropout_rates_ref[i] = (parameters.prior_dropoutrate_mean*parameters.prior_dropoutrate_omega-1+dropoutref_counts[i]) 
                                        / (parameters.prior_dropoutrate_omega-2+nrefalleles[i]);
            dropout_rates_alt[i] = (parameters.prior_dropoutrate_mean*parameters.prior_dropoutrate_omega-1+dropoutalt_counts[i]) 
                                        / (parameters.prior_dropoutrate_omega-2+naltalleles[i]);
            double diff_dropoutrates_lik = dropoutref_counts[i] * std::log(dropout_rates_ref[i]) + (nrefalleles[i] - dropoutref_counts[i]) * std::log(1-dropout_rates_ref[i])
                    + dropoutalt_counts[i] * std::log(dropout_rates_alt[i]) + (naltalleles[i] - dropoutalt_counts[i]) * std::log(1-dropout_rates_alt[i])
                    +(parameters.prior_dropoutrate_mean * parameters.prior_dropoutrate_omega-1) * (std::log(dropout_rates_ref[i]) + std::log(dropout_rates_alt[i])) 
                    +((1.0-parameters.prior_dropoutrate_mean) * parameters.prior_dropoutrate_omega-1) * (std::log(1.0-dropout_rates_ref[i]) + std::log(1.0-dropout_rates_alt[i]));
            double same_dropoutrate_lik = (dropoutref_counts[i]+dropoutalt_counts[i]) * std::log(dropout_rates[i])
                        + (nrefalleles[i] + naltalleles[i] - dropoutref_counts[i] - dropoutalt_counts[i]) * std::log(1-dropout_rates[i])
                        + (parameters.prior_dropoutrate_mean * parameters.prior_dropoutrate_omega-1) * 2*std::log(dropout_rates[i])
                            + ((1.0-parameters.prior_dropoutrate_mean) * parameters.prior_dropoutrate_omega-1) * 2*std::log(1.0-dropout_rates[i]);
            if (same_dropoutrate_lik>diff_dropoutrates_lik-60){
                dropout_rates_ref[i] = dropout_rates[i];
                dropout_rates_alt[i] = dropout_rates[i];
            }
        } 
        else{
            dropout_rates_ref[i] = dropout_rates[i];
            dropout_rates_alt[i] = dropout_rates[i]; 
        }  
        // Set minimum and maximum dropout rate
        if (dropout_rates[i]<0.01) dropout_rates[i]=0.01;
        if (dropout_rates_ref[i]<0.01) dropout_rates_ref[i]=0.01;
        if (dropout_rates_alt[i]<0.01) dropout_rates_alt[i]=0.01;

        if (dropout_rates[i]>0.50) dropout_rates[i]=0.50;
        if (dropout_rates_ref[i]>0.50) dropout_rates_ref[i]=0.50;
        if (dropout_rates_alt[i]>0.50) dropout_rates_alt[i]=0.50;
        
        avg_diff_dropoutrates+= (std::abs(dropout_rates_ref[i] - prev_ref) + std::abs(dropout_rates_alt[i] - prev_alt))  / n_loci;

    }

}

bool Tree::rec_check_max_one_event_per_region_per_lineage(int node, std::vector<int> n_events_in_regions){
    // check that each region is affected, in one lineage, by at most one CNLOH or CNV event.
        for (auto CNV: nodes[node]->get_CNV_events()){
            n_events_in_regions[std::get<0>(CNV)]+=1;
            if (n_events_in_regions[std::get<0>(CNV)]>1) return false;
        }
        for (auto CNLOH: nodes[node]->get_CNLOH_events()){
            n_events_in_regions[std::get<0>(CNLOH)]+=1;
            if (n_events_in_regions[std::get<0>(CNLOH)]>1) return false;
        }
        bool valid=true;
        for (int child: children[node]){
            valid = valid && rec_check_max_one_event_per_region_per_lineage(child,n_events_in_regions);
        }
        return valid;
    }

void Tree::compute_prior_score(){
    // Penalize number of nodes
    log_prior_score=-n_nodes*(4+n_loci)*parameters.node_cost; 
    // Forbid empty nodes
    for (int i=1;i<n_nodes;i++){
        if (nodes[i]->get_number_mutations()==0 && nodes[i]->get_number_CNLOH()==0 && nodes[i]->get_number_CNV()==0) log_prior_score-=100000;
    }

    // Penalize mutations which are not at the root
    log_prior_score+=nodes[0]->get_number_mutations() * parameters.mut_notAtRoot_cost;
    // Particularly penalize mutations with a high pop frequency which are not at the root
    for (int k=1;k<n_nodes;k++){
        for (int mut : nodes[k]->get_mutations()){
            if (data.locus_to_freq[mut]>0.0001) log_prior_score-= data.locus_to_freq[mut] *parameters.mut_notAtRoot_freq_cost;
        }
    }

    // Higher penalty when there are more cells
    double ncells_coef = 0.2 + 1.0*n_cells/8000.0;

    // Penalize CNLOH events (lower penalty at the root, because the germline might be hom alt for SNPs)
    log_prior_score-= ncells_coef*parameters.CNLOH_cost / 10.0 * nodes[0]->get_number_disjoint_CNLOH(); 
    log_prior_score-= ncells_coef*0.1 * nodes[0]->get_number_CNLOH();
    for (int i=1;i<n_nodes;i++){
        log_prior_score-= ncells_coef*parameters.CNLOH_cost * nodes[i]->get_number_disjoint_CNLOH(); 
        log_prior_score-= ncells_coef*0.1 * nodes[i]->get_number_CNLOH(); // slightly higher penalty when there are several CNLOH events on neighbouring regions, to avoid having invalid CNLOH
    }

    // Penalize CNV events
    for (int n=0;n<n_nodes;n++){
        log_prior_score-= ncells_coef*parameters.CNV_cost * nodes[n]->get_number_disjoint_CNV(regions_successor); 
        // Higher penalty for CNVs resulting in LOH, because they have a bigger impact on the likelihood
        log_prior_score-= ncells_coef*parameters.CNV_LOH_cost * nodes[n]->get_number_disjoint_CNV_LOH(regions_successor);  
    }

    // Penalize invalid CNV events
    for (int n=0;n<n_nodes;n++){
        if (!nodes[n]->get_all_CNV_valid()) log_prior_score-= 100000;
    }
    // Cannot have a CNV event at the root.
    if (nodes[0]->get_number_CNV()>0) log_prior_score-= 100000; 
    // One lineage cannot have more than one CNV or CNLOH affecting each region (but it is still possible to have events affecting the same region in parallel branches)
    if (!rec_check_max_one_event_per_region_per_lineage(0,std::vector<int>(n_regions,0))) log_prior_score-=1000000;


    // Dropout rates
    // with probability lambda, both alleles have the same dropout rate. With probability 1-lambda, they have the same dropout rate
    for (int i=0;i<n_loci;i++){
        if (std::abs(dropout_rates_ref[i]-dropout_rates_alt[i])>0.001){
            log_prior_score-=70; 
        }
        // Dropout rates are sampled from a beta distribution.
        log_prior_score+= (parameters.prior_dropoutrate_mean * parameters.prior_dropoutrate_omega-1) 
                                * (std::log(dropout_rates_ref[i]) + std::log(dropout_rates_alt[i])) 
                        + ((1.0-parameters.prior_dropoutrate_mean) * parameters.prior_dropoutrate_omega-1) 
                                * (std::log(1.0-dropout_rates_ref[i]) + std::log(1.0-dropout_rates_alt[i]));
    }
}

void Tree::update_full_score(){
    log_score = log_likelihood + log_prior_score;
}

void Tree::to_dot(std::string filename){
    // Save the tree structure in dot format (for visualization)
    std::vector<std::string> colors{"lightcoral","skyblue3","sandybrown","paleturquoise3","thistle","darkolivegreen3","lightpink","mediumpurple",
                    "darkseagreen3","navajowhite","gold"};
    std::ofstream out_file(filename);

    out_file <<"digraph G{"<<std::endl;
    out_file <<"node [color=dimgray fontsize=24 fontcolor=black fontname=Helvetica penwidth=5];"<<std::endl;
    for (int i=1;i<n_nodes;i++){
        if (parameters.verbose) std::cout<<i<< " is a child of "<<parents[i]<<std::endl;
        out_file<<parents[i]<<" -> "<<i<<" [color=dimgray penwidth=4 weight=2];"<<std::endl;
    }

    for (int i=0;i<n_nodes;i++){
        out_file<<i<<"[label=<"<<nodes[i]->get_label()<<">];"<<std::endl;
        if (parameters.verbose) std::cout<<i<<": "<<nodes[i]->get_label()<<std::endl;
    }

    for (int k=0;k<n_nodes;k++){
        out_file<<k<<" -> "<<k+n_nodes<<" [dir=none style=dashed weight=1 penwidth=5 color="<<colors[k%colors.size()]<<"];"<<std::endl;
    }
    std::vector<int> count_nodes(n_nodes,0);
    int total=0;
    for (int j=0;j<n_cells;j++){
        if (best_attachments[j]>=0 && best_attachments[j]<n_nodes){
            count_nodes[best_attachments[j]]++;
            total++;
        }
    }
    for (int k=0;k<n_nodes;k++){
        double size = std::sqrt(100.0*count_nodes[k]/total) /3.0;
        out_file<<k+n_nodes<<"[label=\""<<count_nodes[k]<<" cells\\n"<<std::round(100.0*count_nodes[k]/total)<<"\\%\""<<" style = filled width="<<size
        <<" height="<<size<<" color="<<colors[k%colors.size()]<<"];"<<std::endl;
    }


    out_file <<"}"<<std::endl;
    out_file.close();
    if (parameters.verbose) std::cout<<"Node probabilities"<<std::endl;
    for (int n=0;n<n_nodes;n++){
        if (parameters.verbose) std::cout<<n<<": "<<node_probabilities[n]<<std::endl;
    }
    if (parameters.verbose) std::cout<<"Dropout rates"<<std::endl;
    for (int i=0;i<n_loci;i++){
        if (parameters.verbose) std::cout<<i<<" ("<<data.locus_to_name[i]<<"): "<<dropout_rates[i]<<" (ref:" <<dropout_rates_ref[i]<<", alt:"<<dropout_rates_alt[i]<<")"<<std::endl;
    }
    
}

void Tree::to_dot_pretty(std::string filename){
    // Save the tree structure in dot format (for visualization)

    std::vector<std::string> colors{"lightcoral","skyblue3","sandybrown","paleturquoise3","thistle","darkolivegreen3","lightpink","mediumpurple",
                    "darkseagreen3","navajowhite","gold"};
    std::ofstream out_file(filename);

    out_file <<"digraph G{"<<std::endl;
    out_file <<"node [color=dimgray fontsize=24 fontcolor=black fontname=Helvetica penwidth=5];"<<std::endl;
    for (int i=1;i<n_nodes;i++){
        if (parameters.verbose) std::cout<<i<< " is a child of "<<parents[i]<<std::endl;
        out_file<<parents[i]<<" -> "<<i<<" [color=dimgray penwidth=4 weight=2];"<<std::endl;
    }

    // Identify mutations at the root which are not affected by a CNV or CNLOH
    std::set<int> excluded_mutations{};
    if (n_nodes>0){
        for (int m: nodes[0]->get_mutations()){
            bool affected_by_event=false;
            for (int n=0;n<n_nodes;n++){
                for (auto CNLOH: nodes[n]->get_CNLOH_events()){
                    if (data.locus_to_region[m]==CNLOH.first) affected_by_event = true;
                }
                for (auto CNV: nodes[n]->get_CNV_events()){
                    if (data.locus_to_region[m]==std::get<0>(CNV)) affected_by_event = true;
                }
            }
            if (!affected_by_event) excluded_mutations.insert(m);
        }
    }

    for (int i=0;i<n_nodes;i++){
        out_file<<i<<"[label=<"<<nodes[i]->get_label_simple(excluded_mutations)<<">];"<<std::endl;
    }

    for (int k=0;k<n_nodes;k++){
        out_file<<k<<" -> "<<k+n_nodes<<" [dir=none style=dashed weight=1 penwidth=5 color="<<colors[k%colors.size()]<<"];"<<std::endl;
    }
    std::vector<int> count_nodes(n_nodes,0);
    int total=0;
    for (int j=0;j<n_cells;j++){
        if (best_attachments[j]>=0 && best_attachments[j]<n_nodes){
             count_nodes[best_attachments[j]]++;
             total++;
        }
    }
    for (int k=0;k<n_nodes;k++){
        double size = std::sqrt(100.0*count_nodes[k]/total) /3.0;
        out_file<<k+n_nodes<<"[label=\""<<count_nodes[k]<<" cells\\n"<<std::round(100.0*count_nodes[k]/total)<<"\\%\""<<" style = filled width="<<size
        <<" height="<<size<<" color="<<colors[k%colors.size()]<<"];"<<std::endl;
    }

    out_file <<"}"<<std::endl;
    out_file.close();
    if (parameters.verbose) std::cout<<"Node probabilities"<<std::endl;
    for (int n=0;n<n_nodes;n++){
        if (parameters.verbose) std::cout<<n<<": "<<node_probabilities[n]<<std::endl;
    }
    if (parameters.verbose) std::cout<<"Dropout rates"<<std::endl;
    for (int i=0;i<n_loci;i++){
        if (parameters.verbose) std::cout<<i<<" ("<<data.locus_to_name[i]<<"): "<<dropout_rates[i]<<" (ref:" <<dropout_rates_ref[i]<<", alt:"<<dropout_rates_alt[i]<<")"<<std::endl;
    }
    
}


Tree::Tree(std::string gv_file, bool use_CNV_arg): //Create tree from a graphviz file
    hastings_ratio(-1.0)     
{
    for (int i=0;i<n_loci;i++){
        dropout_rates.push_back(0.05);
        dropout_rates_ref.push_back(0.05);
        dropout_rates_alt.push_back(0.05);
    }
    cache_scores = new Scores();
    std::ifstream file(gv_file);
    std::string line;
    //skip first 2 lines
    getline (file, line);
    getline (file, line);
    n_nodes=1;
    parents.resize(1);
    parents[0]=-1;
    // Read parents
    bool finished_reading_parents=false;
    while (!finished_reading_parents) {
        getline (file, line);
        if (line.find("->")==std::string::npos) finished_reading_parents=true;
        else{
            int idx=1;
            while (line[idx]!=' ') idx++;
            int parent = stoi(line.substr(0,idx));
            idx++;
            while (line[idx]!=' ') idx++;
            idx++;
            int idx2=idx+1;
            while (line[idx2]!=' '&& line[idx]!=';') idx2++;
            int child = stoi(line.substr(idx,idx2-idx));
            if (child+1 >n_nodes) n_nodes = child+1;
            if (parent+1 > n_nodes) n_nodes = parent+1;
            parents.resize(n_nodes);
            parents[child] = parent;
        }
    }
    // Create nodes 
    for (int i=0;i<n_nodes;i++){
        nodes.push_back(new Node(cache_scores));
    }
    node_probabilities.resize(n_nodes);
    // Read labels (events) and fill the nodes with the events
    bool finished_reading_labels=false;
    while (!finished_reading_labels){
        int idx=1;
        while (idx < line.size() && line[idx]!='[') idx++;
        if (idx>=line.size() || line[idx+1]!='l') finished_reading_labels=true; //empty line
        else{
            int node =  stoi(line.substr(0,idx));
            while (line[idx]!='<') idx++;
            idx++;
            while (line[idx]==' ') idx++;
            bool finished_reading_line=false;
            if (line[idx]=='>') finished_reading_line=true;
            while (!finished_reading_line){
                if (line[idx]=='<') idx+=3; // remove <B>
                if ((line[idx]=='C' && line[idx+2]=='L')){ // CNLOH event
                    idx+=6;
                    int idx2 = idx+1;
                    while (line[idx2]!=':') idx2++;
                    int region = stoi(line.substr(idx,idx2-idx));
                    idx2++;
                    while (line[idx2]!=':') idx2++; // skip region name
                    idx = idx2+1;
                    std::vector<int> lost_alleles{};
                    while (line[idx]!='<' && line[idx]!='b' && line[idx]!='/'){
                        if (line[idx]=='0') lost_alleles.push_back(0);
                        else lost_alleles.push_back(1);
                        idx+=2;
                    }
                    idx--;
                    nodes[node]->add_CNLOH(std::make_pair(region,lost_alleles));
                }
                else if ((line[idx]=='L' && line[idx+2]=='O')){ // CNLOH event
                    idx+=4;
                    int idx2 = idx+1;
                    while (line[idx2]!=':') idx2++;
                    int region = stoi(line.substr(idx,idx2-idx));
                    idx2++;
                    while (line[idx2]!=':') idx2++; // skip region name
                    idx = idx2+1;
                    std::vector<int> lost_alleles{};
                    while (line[idx]!='<' && line[idx]!='b' && line[idx]!='/'){
                        if (line[idx]=='0') lost_alleles.push_back(0);
                        else lost_alleles.push_back(1);
                        idx+=2;
                    }
                    idx--;
                    nodes[node]->add_CNLOH(std::make_pair(region,lost_alleles));
                }
                else if (line[idx]=='C' && line[idx+2]=='V'){ // CNV event
                    int gain_loss=1;
                    if (line[idx+3]=='-') gain_loss=-1;
                    idx+=6;
                    int idx2 = idx+1;
                    while (line[idx2]!=':') idx2++;
                    int region =  stoi(line.substr(idx,idx2-idx));
                    idx = idx2+1;
                    while (line[idx]!=':') idx++;
                    std::vector<int>alleles{};
                    idx++;
                    while (line[idx]!='<' && line[idx]!='b' && line[idx]!='/'){
                        if (line[idx]=='0') alleles.push_back(0);
                        else alleles.push_back(1);
                        idx+=2;
                    }
                    idx--;
                    nodes[node]->add_CNV(std::make_tuple(region,gain_loss,alleles));
                }
                else{ // somatic mutation
                    int idx2= idx+1;
                    while (line[idx2]!=':') idx2++;
                    int locus =  stoi(line.substr(idx,idx2-idx));
                    nodes[node]->add_mutation(locus);
                }
                while (line[idx]!='<') idx++;

                //HTML tags: events separated by <br/>
                if (line[idx+1]=='/') idx+=4; // </B>
                idx+=5;
                if (line[idx]=='>') finished_reading_line=true; 
            }
         } 
        if (!std::getline (file, line)) finished_reading_labels=true;
     }
    // Initialize utils
    cells_attach_loglik.resize(n_cells);
    cells_loglik.resize(n_cells);
    cells_attach_prob.resize(n_cells);
    best_attachments.resize(n_cells);
    use_CNV=false;
    compute_children();
    compute_likelihood(true);
    if(use_CNV_arg && select_regions()){
        use_CNV=true;
        compute_likelihood(true);
    }
    compute_prior_score();
    update_full_score();

    // Close the file
    file.close();
}


void Tree::find_CNV(){
    //1. Find best attachment point of each cell
    int n_attachment_points = n_nodes;
    double doublet_rate = parameters.doublet_rate;
    if (parameters.use_doublets) n_attachment_points = (n_nodes*(n_nodes+3)) / 2; //can attach to a node or to a doublet
    else doublet_rate=0.0;
    double attach_log_lik;
    std::vector<int> best_attachments(n_cells,-1);
    double Z;
    

    std::vector<double> attach_prob{}; // Posterior probability of the attachments of a cell
    attach_prob.resize(n_attachment_points);
    std::vector<std::vector<double>> cells_attach_prob{};
    for (int j=0; j<n_cells; j++){ 
        attach_prob.clear();
        double best_attach_score=-DBL_MAX;
        int best_attachment=-1;
        for (int k=0; k<n_nodes;k++){
            attach_log_lik = nodes[k]->attachment_scores[j] + std::log(node_probabilities[k]) +std::log(1-doublet_rate);
            if (attach_log_lik>best_attach_score) {
                best_attach_score=attach_log_lik;
                best_attachment=k;
            }
        }
        if (parameters.use_doublets){
            int idx=0;
            for (int k=0;k<n_nodes;k++){
                for (int l=k;l<n_nodes;l++){
                    attach_log_lik = doublets[idx]->attachment_scores[j] + std::log(node_probabilities[k]) 
                    + std::log(node_probabilities[l]) +std::log(doublet_rate);
                    if (k!=l)  attach_log_lik+= std::log(2); // the doublet (l,k) has the same probability as (k,l)
                    idx++;
                    if (attach_log_lik>best_attach_score) best_attachment=-1;
                }
            }
        }
        best_attachments[j] = best_attachment;
    }

    //2. Compute average region probability in each node, and find if there is a difference between some nodes
    for (int i=0;i<n_regions;i++){
        std::vector<std::vector<double>> nodes_regionprobs{};
        std::vector<double> nodes_averages(n_nodes,0.0);
        nodes_regionprobs.resize(n_nodes);
        for (int j=0;j<n_cells;j++){
            if (best_attachments[j]>=0){
                nodes_regionprobs[best_attachments[j]].push_back(1.0*cells[j].region_counts[i] / cells[j].total_counts);
            }
        }
        double min_prop=+DBL_MAX;
        double max_prop=-DBL_MAX;
        int min_node=-1;
        int max_node=-1;
        for (int k=0;k<n_nodes;k++){
            if (nodes_regionprobs[k].size()>8){
                double prop = 0;
                for (double d: nodes_regionprobs[k]) prop+= d/nodes_regionprobs[k].size();
                nodes_averages[k] = prop;
                if (prop>max_prop){
                    max_prop = prop;
                    max_node = k;
                }
                if (prop<min_prop){
                    min_prop = prop;
                    min_node = k;
                }
            }
        }
        //if (max_prop>min_prop*1.1 && max_prop>0.1/n_regions){
        if (true){
            std::cout<<"candidate CNV in region "<<i<<" ("<<data.region_to_name[i]
            <<"). Node "<<min_node<<": avg "<<min_prop<<", node "<<max_node<<": avg "<<max_prop <<" (ratio "<<max_prop/min_prop<<")"<<std::endl;
            for (int k=0;k<n_nodes;k++){
                if (nodes_regionprobs[k].size()>8) std::cout<<"node "<<k<<": size "<<nodes_regionprobs[k].size()<<", avg "<<nodes_averages[k]<<std::endl;
            }
        }
    }
}

bool Tree::select_regions(int index){
    // Find regions which might contain a CNV event. Return true if it was possible to estimate node regions (if the node contains enough cells), false otherwise.
    candidate_regions.clear();
    region_probabilities.resize(n_regions);
    //  Compute average region probability in each node, and find if there is a difference between some nodes
    for (int k=0;k<n_regions;k++){
        if (!data.region_is_reliable[k]) continue;
        std::vector<std::vector<double>> nodes_regionprobs{};
        std::vector<double> nodes_averages(n_nodes,0.0);
        nodes_regionprobs.resize(n_nodes);
        for (int j=0;j<n_cells;j++){
            if (best_attachments[j]<n_nodes && best_attachments[j]>=0){
                nodes_regionprobs[best_attachments[j]].push_back(1.0*cells[j].region_counts[k] / cells[j].total_counts);
            }
        }
        if (nodes_regionprobs[0].size()<std::max(40.0,0.015*n_cells)){
            if (index >=0) std::cout<<"Chain "<<std::to_string(index)<<": ";
            std::cout<<"In the tree inferred without CNVs, there were not enough cells attached to the root to estimate the region weights, so COMPASS could not attempt to find CNVs."<<std::endl;
            return false; // not enough cells attached to the root: cannot find CNVs
        }
        double rootprob=0;
        for (double prob: nodes_regionprobs[0]){
            rootprob+= prob / nodes_regionprobs[0].size();
            if (data.region_to_chromosome[k]!="X" || data.sex=="female") region_probabilities[k] = rootprob;
            else region_probabilities[k] = rootprob*2; // region probability is defined for diploid
        }

        // Select regions whose probability is different between the root and another node
        for (int n=1;n<n_nodes;n++){
            if (nodes_regionprobs[n].size()>=std::max(40.0,0.03*n_cells)){
                double prob = 0;
                for (double d: nodes_regionprobs[n]) prob+= d/nodes_regionprobs[n].size();
                if ((prob>rootprob*1.275 || rootprob>prob*1.35) && rootprob>0.05/n_regions){
                    candidate_regions.push_back(k);
                    break;
                }
            }
        }
    }
    if (candidate_regions.size()==0){
        if (index >=0) std::cout<<"Chain "<<std::to_string(index)<<": ";
        std::cout<<"In the tree inferred without CNVs, there were no regions whose coverage in one node were different from the root, so COMPASS could not identify any CNVs."<<std::endl;
        return false;
    }


    regions_successor = std::vector<int>(n_regions,-1);
    for (int k=0;k<candidate_regions.size()-1;k++){
        int k1 = candidate_regions[k];
        int k2 = candidate_regions[k+1];
        if (data.region_to_chromosome.size()>0 && data.region_to_chromosome[k1] == data.region_to_chromosome[k2]){
            regions_successor[k1]=k2;
        }
    }
    if (true || parameters.verbose){
        if (index >=0) std::cout<<"Chain "<<std::to_string(index)<<": ";
        std::cout<<"After the first phase (inferring the best tree without CNVs), COMPASS identified the following candidate regions which might contain CNVs: ";
        for (int i: candidate_regions){
            std::cout<<data.region_to_name[i]<<",";
        }
        std::cout<<std::endl;
    }
    
    return true;
}


void Tree::add_node(int parent){
    // Add a new node below an existing one, and randomly reassign the children of the parent node.
    n_nodes++;
    nodes.push_back(new Node(cache_scores));
    parents.push_back(parent);
    node_probabilities.resize(n_nodes);
    for (int child: children[parent]){
        if (std::rand()%2==0) parents[child] = n_nodes-1;
    }
    compute_children();
}

void Tree::delete_node(int node){
    for (int child: children[node]){
        parents[child] = parents[node];
    }
    delete nodes[node];
    // Change the index of the last node.
    if (node!=n_nodes-1){
        nodes[node] = nodes[n_nodes-1];
        for (int n=0;n<n_nodes;n++){
            if (parents[n]==n_nodes-1) parents[n] = node;
        }
        parents[node] = parents[n_nodes-1];
    }
    n_nodes--;
    nodes.resize(n_nodes);
    parents.resize(n_nodes);
    node_probabilities.resize(n_nodes);
    compute_children();
}


void Tree::prune_reattach(){
    // Prune a subtree, and reattach it somewhere else in the tree

    if (n_nodes<=1){ // Need at least 2 nodes to be able to prune and reattach.
        hastings_ratio=0.0;
        return;
    }

    // Randomly select one node (apart from the root) to prune and reattach
    int pruned_node_id = ( std::rand() % (n_nodes-1) ) + 1;

    // Find all nodes that are not a descendant of the pruned node
    // Perform a DFT from the pruned node to identify all its descendants
    std::vector<bool> is_descendant(n_nodes,false);
    std::stack<int> stk;
    stk.push(pruned_node_id);
    while (!stk.empty()) {
        int top = stk.top();
        stk.pop();
        for (int child:children[top]) {
            stk.push(child);
        }
        is_descendant[top] = true;
    }
    std::vector<int> candidate_attachment_points;
    for (int i=0;i<n_nodes;i++){
        if (!is_descendant[i])
            candidate_attachment_points.push_back(i);
    }

    // Sample one attachment point among all nodes that are not descendant of the pruned node
    int attachment_point_id = candidate_attachment_points[std::rand()%candidate_attachment_points.size()];

    // Update the tree accordingly, and recompute the genotypes
    parents[pruned_node_id] = attachment_point_id;
    compute_children(); //update the children list according to the new parent vector

    hastings_ratio=1.0;
}

void Tree::swap_node_labels(){
    // Select the 2 nodes and copy the nodes below them
    if (n_nodes<=1){
        hastings_ratio=0.0;
        return;
    }
    int node1 = std::rand()%n_nodes;
    int node2 = std::rand()%(n_nodes-1);
    if (node2==node1) node2 = n_nodes-1;
    // Exchange the nodes
    Node* temp = nodes[node1];
    nodes[node1] = nodes[node2];
    nodes[node2] = temp;

    hastings_ratio=1.0;
}

/*void Tree::add_remove_mutation(){
    // Find nodes with mutations
    std::vector<int> nodes_with_mut{};
    for (int i=0;i<n_nodes;i++){
        if (nodes[i]->get_number_mutations()>0) nodes_with_mut.push_back(i);
    }
    // Find mutations missing from the tree
    std::set<int> mutations_missing{};
    for (int i=0;i<n_loci;i++) mutations_missing.insert(i);
    for (int i=0;i<n_nodes;i++){
        for (int k: nodes[i]->get_mutations()) mutations_missing.erase(k);
    }
    //for (int k: mutations_missing) std::cout<<"mut "<<k<<" missing"<<std::endl;
    double add_probability = 1.0*mutations_missing.size()/n_loci;
    if (1.0*std::rand()/RAND_MAX<add_probability){
        // Add mutation
        //std::cout<<"Add mutation"<<std::endl;
        // Select mutation to add among the missing ones
        int mut = *std::next(mutations_missing.begin(), std::rand()%mutations_missing.size());
        // Select node where to add the mutation
        int node = std::rand()%n_nodes;
        nodes[node]->add_mutation(mut);

        // To reverse the move, we need to
        // 1. select remove
        // 2. select this node
        // 3. select this mutation
        
        int n_nodes_with_mut= nodes_with_mut.size();
        // if we added a mutation to a node that had no mutation, there is now one more node with mutations
        if (nodes[node]->get_number_mutations()==1) n_nodes_with_mut++;
        hastings_ratio = (1.0-add_probability)/add_probability *n_nodes / n_nodes_with_mut / nodes[node]->get_number_mutations() * mutations_missing.size();
    }
    else{
        // Remove mutation
        //std::cout<<"Remove mutation"<<std::endl;
        // Select one node with at least a mutation
        int node = nodes_with_mut[std::rand()%nodes_with_mut.size()];
        int mut = nodes[node]->remove_mutation();
        double new_add_prob = (1.0+mutations_missing.size()) / n_nodes;;

        hastings_ratio = new_add_prob / (1.0-add_probability) * nodes_with_mut.size() / n_nodes 
                        * (nodes[node]->get_number_mutations()+1.0) / (mutations_missing.size()+1);
    }

}*/

void Tree::move_mutation(){

    // sample the source node (from which the event will be moved) among the nodes which have at least one mutation.
    std::vector<int> nodes_with_event{};
    for (int i=0;i<n_nodes;i++){
        if (nodes[i]->get_number_mutations()>0) nodes_with_event.push_back(i);
    }
    int initial_nb_nodes_with_events = nodes_with_event.size();
    int new_nb_nodes_with_events = initial_nb_nodes_with_events;
    int source_node = nodes_with_event[std::rand() % nodes_with_event.size()];
    hastings_ratio=1.0 * initial_nb_nodes_with_events * nodes[source_node]->get_number_mutations();

    // Select the destination node (either an existing or a new node)
    int destination_node;
    double new_node_prob=0.2;
    if (n_nodes<=1) new_node_prob = 1.0;
    if (1.0*std::rand()/RAND_MAX<new_node_prob){
        // Move the event to a new node.
        int parent = std::rand()%n_nodes; // Select the parent of the new node
        hastings_ratio*= 1/new_node_prob * n_nodes * std::pow(2,children[parent].size());
        add_node(parent);
        destination_node = n_nodes-1;
    }
    else{
        // Move the event to an existing node
        destination_node = std::rand()%(n_nodes-1);
        if (destination_node==source_node) destination_node = n_nodes-1;
        hastings_ratio*= 1/(1.0-new_node_prob) * (n_nodes-1);
        
    }

    // Remove the mutation from the source node and add it to the destination node
    int mutation = nodes[source_node]->remove_random_mutation();
    nodes[destination_node]->add_mutation(mutation);
    if (nodes[destination_node]->get_number_mutations()==1) new_nb_nodes_with_events++;

    // If the source node is now empty, remove it.
    if (source_node !=0 && nodes[source_node]->is_empty()){
        // Delete the source node
        // To reverse the move, we will need to select add new node, select the right parent for the new destination node, and reassign the children correctly.
        hastings_ratio*= new_node_prob / n_nodes/ std::pow(2,children[source_node].size());
        delete_node(source_node);
        new_nb_nodes_with_events--;
    }
    else{
        // The source node still exists. To reverse the move, we need to select not adding a new node, and select the right destination node.
        hastings_ratio*= (1.0-new_node_prob) / (n_nodes-1);
        if (nodes[source_node]->get_number_mutations()==0) new_nb_nodes_with_events--;
    }

    // To reverse the move, we must select the right source node, the right event in this node
    hastings_ratio*= 1.0/new_nb_nodes_with_events / nodes[destination_node]->get_number_mutations();
}

void Tree::split_merge_node(){
    double default_merge_probability = 0.20;
    double merge_probability = default_merge_probability;
    if (n_nodes<=1) merge_probability =0.0; //cannot merge nodes if there is only one node.

    if ( (1.0*std::rand())/RAND_MAX <= merge_probability){ // merge nodes
        // Select one node which is not the root
        int node1 = std::rand()%(n_nodes-1) + 1;
        int node2 = parents[node1];
        // Move the events of node 1 into node 2 
        while (nodes[node1]->get_number_mutations()>0){
            nodes[node2]->add_mutation(nodes[node1]->remove_random_mutation());
        }
        while (nodes[node1]->get_number_CNLOH()>0){
            std::pair<int,std::vector<int>> CNLOH = nodes[node1]->remove_random_CNLOH();
            nodes[node2]->add_CNLOH(CNLOH);
        }
        while (nodes[node1]->get_number_CNV()>0){
            std::tuple<int,int,std::vector<int>> CNV = nodes[node1]->remove_random_CNV();
            nodes[node2]->add_CNV(CNV);
        }

        delete_node(node1);
        if ( node1!=n_nodes-1 && node2==n_nodes-1) node2=node1; // the index of node 2 might have changed following the deletion.

        // Hastings ratio
        // For merge, there is one possibility to merge the events and set the new parents
        // For split, there are 2**n_events possibilities to split the events between the 2 nodes
        // and 2**n_children possibilities to set the parent of the children of the split node
        if (n_nodes==1) hastings_ratio = 1.0/merge_probability;
        else hastings_ratio = (1.0-merge_probability)/merge_probability;
        
        hastings_ratio *= 1.0/std::pow(2.0,nodes[node2]->get_number_mutations() + nodes[node2]->get_number_CNLOH() + nodes[node2]->get_number_CNV());
        hastings_ratio *= 1.0/ std::pow(2.0,children[node2].size()); 
    }
    else{ // Split node
        int node = std::rand()%n_nodes;
        add_node(node); // The new node is a child of the original node, and the children of the original node are randomly distributed between the 2 nodes.

        //Randomly split the events between the two nodes (choosing a random number of events to move)
        if (nodes[node]->get_number_mutations()>0){
            int n_events_moved = std::rand()%(nodes[node]->get_number_mutations()+1) ; // we can move from 0 to all of the mutations
            for (int i=0;i<n_events_moved;i++){
                nodes[n_nodes-1]->add_mutation(nodes[node]->remove_random_mutation());
            }
        }
        if (nodes[node]->get_number_CNLOH()>0){
            int n_events_moved = std::rand()%(nodes[node]->get_number_CNLOH()+1) ;
            for (int i=0;i<n_events_moved;i++){
                std::pair<int,std::vector<int>> CNLOH = nodes[node]->remove_random_CNLOH();
                nodes[n_nodes-1]->add_CNLOH(CNLOH);
            }
        }
        if (nodes[node]->get_number_CNV()>0){
            //TODO: maybe make it more likely to move adjacent CNVs together ?
            int n_events_moved = std::rand()%(nodes[node]->get_number_CNV()+1) ;
            for (int i=0;i<n_events_moved;i++){
                std::tuple<int,int,std::vector<int>> CNV = nodes[node]->remove_random_CNV();
                nodes[n_nodes-1]->add_CNV(CNV);
            }
        }

        // Hastings ratio
        int n_events = nodes[node]->get_number_mutations()+nodes[n_nodes-1]->get_number_mutations()
                        + nodes[node]->get_number_CNLOH()+nodes[n_nodes-1]->get_number_CNLOH()
                        + nodes[node]->get_number_CNV()+nodes[n_nodes-1]->get_number_CNV();
        if (n_nodes==2) hastings_ratio = default_merge_probability;
        else hastings_ratio = merge_probability / (1.0-merge_probability);
        
        hastings_ratio *= std::pow(2.0,n_events);
        hastings_ratio *= std::pow(2.0,children[node].size()+children[n_nodes-1].size()); 
    }
}



void Tree::add_remove_CNLOH(){
    double default_add_probability=0.80;
    double add_probability = default_add_probability;
    // Find nodes which have loh events
    
    std::vector<int> nodes_with_events{};
    for (int i=0;i<n_nodes;i++){
        if (nodes[i]->get_number_CNLOH()>0) nodes_with_events.push_back(i);
    }
    
    if (nodes_with_events.size()==0) add_probability=1.0; //cannot remove a cnloh event if none exists

    std::vector<int> regions_with_loci{};
    for (int k=0;k<data.region_to_loci.size();k++){
        if (data.region_to_loci[k].size()>0) regions_with_loci.push_back(k);
    }

    if ( (1.0*std::rand())/RAND_MAX <= add_probability){ // add CNLOH event
        // Select node, locus and which alleles are lost
        int node = std::rand()%(2*n_nodes);
        int region = regions_with_loci[std::rand()%regions_with_loci.size()];
        std::vector<int> lost_alleles;
        for (int i=0;i<data.region_to_loci[region].size();i++) lost_alleles.push_back(std::rand()%2);
        if (node<n_nodes){
            // Add CNLOH event to an existing node
            hastings_ratio=1.0;
        }
        else{
            // Add CNLOH event to a new node
            int parent = node - n_nodes;
            node = n_nodes;
            add_node(parent);
            hastings_ratio = std::pow(2,children[parent].size());
        }
        nodes[node]->add_CNLOH(std::make_pair(region,lost_alleles));
        //There are n_nodes possibilities for where to place the CNLOH, n_region_with_loci possibilities for the region and 2**(n_loci in region) possibilities for the alleles
        // To reverse the move, we need to select remove, select the same node, and select the right CNLOH event
        int n_nodes_with_events = nodes_with_events.size();
        if (nodes[node]->get_number_CNLOH()==1) n_nodes_with_events++;
        hastings_ratio *= (1.0-default_add_probability) /n_nodes_with_events  / nodes[node]->get_number_CNLOH() 
                            / add_probability * 2 * n_nodes * regions_with_loci.size() * std::pow(2,data.region_to_loci[region].size());
    }
    else{ // remove CNLOH event
        // Select a node which has a CNLOH event
        int node = nodes_with_events[std::rand()%nodes_with_events.size()];
        std::pair<int,std::vector<int>> CNLOH = nodes[node]->remove_random_CNLOH(); // randomly remove one loh event
        hastings_ratio=1.0;
        if (node!=0 && nodes[node]->is_empty()){
            // The node is now empty--> delete it.
            // Set the parent of the children of the deleted node.
            int parent = parents[node];
            if (node!=n_nodes-1 && parent==n_nodes-1) parent = node;
            delete_node(node);
            hastings_ratio = 1.0/std::pow(2,children[parent].size());
        }
        hastings_ratio *= add_probability / n_nodes / regions_with_loci.size() / std::pow(2,data.region_to_loci[CNLOH.first].size())
                        / (1.0-add_probability) * nodes_with_events.size() * (nodes[node]->get_number_CNLOH()+1);
    }
}

void Tree::move_CNLOH(){
    // sample the source node (from which the event will be moved) among the nodes which have at least one CNV event.
    std::vector<int> nodes_with_event{};
    for (int i=0;i<n_nodes;i++){
        if (nodes[i]->get_number_CNLOH()>0) nodes_with_event.push_back(i);
    }
    int initial_nb_nodes_with_events = nodes_with_event.size();
    int new_nb_nodes_with_events = initial_nb_nodes_with_events;
    if (initial_nb_nodes_with_events==0){
        hastings_ratio=0.0;
        return;
    }
    double new_node_prob=0.20;
    int source_node = nodes_with_event[std::rand() % nodes_with_event.size()];
    hastings_ratio=1.0 * initial_nb_nodes_with_events * nodes[source_node]->get_number_CNLOH();
    std::pair<int,std::vector<int>> CNLOH = nodes[source_node]->remove_random_CNLOH();
    if (source_node !=0 && nodes[source_node]->is_empty()){
        // Delete the source node
        // To reverse the move, we will need to select add new node, select the right parent for the new destination node, and reassign the children correctly.
        hastings_ratio*= new_node_prob / n_nodes/ std::pow(2,children[source_node].size());
        delete_node(source_node);
        new_nb_nodes_with_events--;
    }
    else{
        // The source node still exists. To reverse the move, we need to select not adding a new node, and select the right destination node.
        hastings_ratio*= (1.0-new_node_prob) / (n_nodes-1);
        if (nodes[source_node]->get_number_CNLOH()==0) new_nb_nodes_with_events--;
    }
    
    int destination_node;
    if (n_nodes<=1) new_node_prob = 1.0;
    if (1.0*std::rand()/RAND_MAX<new_node_prob){
        // Move the event to a new node.
        int parent = std::rand()%n_nodes; // Select the parent of the new node
        hastings_ratio*= 1/new_node_prob * n_nodes * std::pow(2,children[parent].size());
        add_node(parent);
        destination_node = n_nodes-1;
    }
    else{
        // Move the event to an existing node
        destination_node = std::rand()%(n_nodes-1);
        if (destination_node==source_node) destination_node = n_nodes-1;
        hastings_ratio*= 1/(1.0-new_node_prob) * (n_nodes-1);
    }
    // Add the event, potentially altering the sign
    nodes[destination_node]->add_CNLOH(CNLOH);
    if (nodes[destination_node]->get_number_CNLOH()==1) new_nb_nodes_with_events++;

    // To reverse the move, we must select the right source node, the right event in this node
    hastings_ratio*= 1.0/new_nb_nodes_with_events / nodes[destination_node]->get_number_CNLOH();
}


void Tree::add_remove_CNV(){
    double default_add_probability=0.70;
    double add_probability = default_add_probability;
    
    std::vector<int> nodes_with_events{};
    for (int i=0;i<n_nodes;i++){
        if (nodes[i]->get_number_CNV()>0) nodes_with_events.push_back(i);
    }
    if (candidate_regions.size()==0){
        hastings_ratio=0.0;
        return;
    }
    if (nodes_with_events.size()==0)  add_probability=1.0; //cannot remove a CNV event if none exists
    

    if ( (1.0*std::rand())/RAND_MAX <= add_probability){ // add CNV event
        // Select node, region, gain/loss and alleles
        int node = std::rand()%(2*n_nodes); // can add the mutation to an existing node, or to a new node below an existing node
        int region = candidate_regions[std::rand()%candidate_regions.size()];
        int gain_loss = 1 - 2*(std::rand()%2);
        int parent;
        if (node<n_nodes){ // Add event to an existing node
            parent = node;
            hastings_ratio=1.0;
        }
        else{ // Create a new node
            parent = node - n_nodes;
            node = n_nodes;
            add_node(parent);
            hastings_ratio = std::pow(2,children[parent].size()); 
        }
        std::vector<int> alleles{};
        for (int i=0;i<data.region_to_loci[region].size();i++){
            int allele =std::rand()%2; //0: CNV affects ref allele; 1: CNV affects alt allele
            // Can only gain/lose lose one allele if we had at least one copy of it !
            if (allele==0 && nodes[parent]->get_n_ref_allele(data.region_to_loci[region][i])==0) allele=1;
            if (allele==1 && nodes[parent]->get_n_alt_allele(data.region_to_loci[region][i])==0) allele=0;
            alleles.push_back(allele);
        }
        nodes[node]->add_CNV(std::make_tuple(region,gain_loss,alleles));
        // There are 2*n_nodes possibilities for where to place the CNV, n_regions possibilities for the regions,2 possibilities for gain/loss,
        // 2**n_alleles possibilities for the alleles and in case a new node was created, 2**(n_children of the parent) possibilities to assign the children of the parent
        // To reverse the move, we need to select remove, select the same node, and select the right CNV event
        int n_nodes_with_events = nodes_with_events.size();
        if (nodes[node]->get_number_CNV()==1) n_nodes_with_events++;
        hastings_ratio *= (1.0-default_add_probability) /n_nodes_with_events  / nodes[node]->get_number_CNV() 
                            / add_probability * 2*n_nodes * candidate_regions.size() * 2 * std::pow(2,alleles.size());
    }
    else{ // remove CNV event
        // Select a node which has a CNV event
        int node = nodes_with_events[std::rand()%nodes_with_events.size()];
        auto CNV = nodes[node]->remove_random_CNV(); // randomly remove one CNV event
        int n_alleles = std::get<2>(CNV).size();
        hastings_ratio=1.0;
        if (node!=0 && nodes[node]->is_empty()){
            // The node is now empty--> delete it.
            // Set the parent of the children of the deleted node.
            int parent = parents[node];
            if (node!=n_nodes-1 && parent==n_nodes-1) parent = node;
            delete_node(node);
            hastings_ratio = 1.0/std::pow(2,children[parent].size());
        }
        hastings_ratio *= add_probability / 2.0 / n_nodes / candidate_regions.size() /2.0 / std::pow(2,n_alleles)
                            / (1.0-add_probability) * nodes_with_events.size() * (nodes[node]->get_number_CNV()+1);
    }
}


void Tree::move_CNV(){
    // sample the source node (from which the event will be moved) among the nodes which have at least one CNV event.
    std::vector<int> nodes_with_event{};
    for (int i=0;i<n_nodes;i++){
        if (nodes[i]->get_number_CNV()>0) nodes_with_event.push_back(i);
    }
    int initial_nb_nodes_with_events = nodes_with_event.size();
    int new_nb_nodes_with_events = initial_nb_nodes_with_events;
    if (initial_nb_nodes_with_events==0){
        hastings_ratio=0.0;
        return;
    }
    
    double new_node_prob=0.20;
    int source_node = nodes_with_event[std::rand() % nodes_with_event.size()];
    hastings_ratio=1.0 * initial_nb_nodes_with_events * nodes[source_node]->get_number_CNV();
    std::tuple<int,int,std::vector<int>> CNV = nodes[source_node]->remove_random_CNV();
    if (source_node !=0 && nodes[source_node]->is_empty()){
        // Delete the source node
        // To reverse the move, we will need to select add new node, select the right parent for the new destination node, and reassign the children correctly.
        hastings_ratio*= new_node_prob / n_nodes/ std::pow(2,children[source_node].size());
        delete_node(source_node);
        new_nb_nodes_with_events--;
    }
    else{
        // The source node still exists. To reverse the move, we need to select not adding a new node, and select the right destination node.
        hastings_ratio*= (1-new_node_prob) / (n_nodes-1);
        if (nodes[source_node]->get_number_CNV()==0) new_nb_nodes_with_events--;
    }
    
    int destination_node;
    if (n_nodes<=1) new_node_prob = 1.0;
    if (1.0*std::rand()/RAND_MAX<new_node_prob){
        // Move the event to a new node.
        int parent = std::rand()%n_nodes; // Select the parent of the new node
        hastings_ratio*= 1.0/new_node_prob * n_nodes * std::pow(2,children[parent].size());
        add_node(parent);
        destination_node = n_nodes-1;
    }
    else{
        // Move the event to an existing node
        destination_node = std::rand()%(n_nodes-1);
        if (destination_node==source_node) destination_node = n_nodes-1;
        hastings_ratio*= 1.0/(1.0-new_node_prob) * (n_nodes-1);
    }
    // Add the event
    nodes[destination_node]->add_CNV(CNV);
    if (nodes[destination_node]->get_number_CNV()==1) new_nb_nodes_with_events++;

    // To reverse the move, we must select the right source node, the right event in this node
    hastings_ratio*= 1.0/new_nb_nodes_with_events / nodes[destination_node]->get_number_CNV();
}

void Tree::merge_or_duplicate_CNV(){
    // If there are two indentical CNVs in the tree, can merge them into one CNV at their most recent common ancestor
    // If one node containing a CNV has multiple children, can duplicate this CNV and place both copies in parallel branches.

    // Find duplicate CNVs
    std::multiset<std::tuple<int,int,std::vector<int>>> CNVs_in_tree;
    for (int n=0;n<n_nodes;n++){
        for (auto CNV: nodes[n]->get_CNV_events()) CNVs_in_tree.insert(CNV);
    }
    std::vector<std::tuple<int,int,std::vector<int>>> duplicate_CNVs{};
    for (auto CNV: CNVs_in_tree){
        if (CNVs_in_tree.count(CNV)>1){
            // Only add the CNV if it was not already in the vector
            bool new_CNV=true;
            for (auto CNV_already_in_vector: duplicate_CNVs){
                new_CNV = new_CNV && (CNV!=CNV_already_in_vector);
            }
            if (new_CNV) duplicate_CNVs.push_back(CNV);
        }
    }

    // Find nodes containing a CNV event and having multiple children
    std::vector<int> nodes_with_CNV_and_multiple_children{};
    for (int n=0;n<n_nodes;n++){
        if (nodes[n]->get_number_CNV()>0 && children[n].size()>1) nodes_with_CNV_and_multiple_children.push_back(n);
    }

    if (duplicate_CNVs.size()==0 && nodes_with_CNV_and_multiple_children.size()==0){
        // No CNV can be merged or duplicated
        hastings_ratio=0.0;
        return;
    } 


    int event_index = std::rand()%(duplicate_CNVs.size() + nodes_with_CNV_and_multiple_children.size());
    if (event_index < duplicate_CNVs.size()){
        // Merge 2 CNVs
        auto CNV = duplicate_CNVs[event_index];
        std::vector<int> nodes_containing_CNV{};
        for (int n=0;n<n_nodes;n++){
            for (auto CNV_node: nodes[n]->get_CNV_events()){
                if (CNV==CNV_node) nodes_containing_CNV.push_back(n);
            }
        }
        int index1 = std::rand()%nodes_containing_CNV.size();
        int index2 = std::rand()%(nodes_containing_CNV.size()-1);
        if (index1==index2) index2=nodes_containing_CNV.size()-1;
        int node1 = nodes_containing_CNV[index1];
        int node2 = nodes_containing_CNV[index2];
        // Find their most recent common ancestor
        int ancestor = node1;
        while (!is_ancestor(ancestor,node2)) ancestor = parents[ancestor];
        nodes[node1]->remove_CNV(CNV);
        nodes[node2]->remove_CNV(CNV);
        int new_n_nodes_with_CNV_and_multiple_children = nodes_with_CNV_and_multiple_children.size();
        if (std::rand()%2==0){
            add_node(ancestor);
            nodes[n_nodes-1]->add_CNV(CNV);
            new_n_nodes_with_CNV_and_multiple_children+=1;
        }
        else{
            nodes[ancestor]->add_CNV(CNV);
            if (nodes[ancestor]->get_number_CNV()==1) new_n_nodes_with_CNV_and_multiple_children+=1;
        }

        // Hastings ratio
        // In order to merge, we had to:
        //  * Select merge (probability: duplicate_CNVs.size() / (duplicate_CNVs.size() + nodes_with_CNV_and_multiple_children.size()))
        //  * Select the CNV among duplicate_CNVs.size() possibilities
        //  * Select the 2 nodes containing this CNV: nodes_containing_CNV.size() choose 2 possibilities
        //  * Select whether to move the CNV to their MRCA or create a new node below it (2 possibilities)
        // In order to reverse this move, we would have to:
        //  * Select duplicate (new_n_nodes_with_CNV_and_multiple_children / (duplicate_CNVs.size()-1 + new_n_nodes_with_CNV_and_multiple_children))
        //       (new_n_nodes_with_CNV_and_multiple_children because after having performed the move, there might be one more node with CNV and multiple children)
        //  * Select the right node among nodes_with_CNV_and_multiple_children_new.size() possibilities
        //  * Select the right CNV among nodes[ancestor]->get_number_CNV() possibilities
        //  * Select the right 2 children: children[ancestor].size() possibilities
        //  * For each child, select the right descendant to put the event.

        int n_descendants1=0; // Number of descendants in the subtree below "ancestor" containing node1
        int node=node1;
        while (node!=ancestor){
            node = parents[node];
            n_descendants1+=1;
        }
        std::stack<int> stk;
        stk.push(node1);
        while (!stk.empty()) {
            int top = stk.top();
            stk.pop();
            for (int child: children[top]) stk.push(child);
            if (top!=node1) n_descendants1++;
        }
        int n_descendants2=0;
        node=node2;
        while (node!=ancestor){
            node = parents[node];
            n_descendants2+=1;
        }
        stk.push(node2);
        while (!stk.empty()) {
            int top = stk.top();
            stk.pop();
            for (int child: children[top]) stk.push(child);
            if (top!=node2) n_descendants2++;
        }
        // Select merge
        hastings_ratio= 1.0 / duplicate_CNVs.size() * (duplicate_CNVs.size() + nodes_with_CNV_and_multiple_children.size());
        // Select which CNVs to merge, the nodes containing that CNV and whether to put the CNV to the MRCA or a new node below it.
        hastings_ratio*= duplicate_CNVs.size() * std::exp(cache_scores->log_n_choose_k(nodes_containing_CNV.size(),2) *2) *2.0;
        // Reverse: Select duplicate
        hastings_ratio*= 1.0* new_n_nodes_with_CNV_and_multiple_children / (duplicate_CNVs.size()-1 + new_n_nodes_with_CNV_and_multiple_children);
        // Select the node, the CNV and the descendants that will get a copy of the CNV
        hastings_ratio*= 1.0 / new_n_nodes_with_CNV_and_multiple_children / nodes[ancestor]->get_number_CNV() / children[ancestor].size() / n_descendants1 / n_descendants2;
    }
    else{
        // Duplicate 1 CNV

        //Select the node and the CNV event to duplicate
        int node_ancestor = nodes_with_CNV_and_multiple_children[std::rand()%nodes_with_CNV_and_multiple_children.size()];
        std::vector<std::tuple<int,int,std::vector<int>>> CNVs_node{};
        for (auto CNV: nodes[node_ancestor]->get_CNV_events()) CNVs_node.push_back(CNV);
        auto CNV = CNVs_node[std::rand()%CNVs_node.size()];

        // Select the 2 nodes where to add the CNV
        int child1=children[node_ancestor][std::rand()%children[node_ancestor].size()];
        int child2=children[node_ancestor][std::rand()%children[node_ancestor].size()];
        while (child1==child2) child2=children[node_ancestor][std::rand()%children[node_ancestor].size()];

        std::vector<int> descendants1{};
        std::stack<int> stk;
        stk.push(child1);
        while (!stk.empty()) {
            int top = stk.top();
            stk.pop();
            for (int child: children[top]) stk.push(child);
            descendants1.push_back(top);
        }
        std::vector<int> descendants2{};
        stk.push(child2);
        while (!stk.empty()) {
            int top = stk.top();
            stk.pop();
            for (int child: children[top]) stk.push(child);
            descendants2.push_back(top);
        }
        int node1 = descendants1[std::rand()%descendants1.size()];
        int node2 = descendants2[std::rand()%descendants2.size()];

        nodes[node_ancestor]->remove_CNV(CNV);
        if (node_ancestor>0 && nodes[node_ancestor]->is_empty()) delete_node(node_ancestor);
        nodes[node1]->add_CNV(CNV);
        nodes[node2]->add_CNV(CNV);

        int new_n_nodes_with_CNV_and_multiple_children = nodes_with_CNV_and_multiple_children.size();
        if (nodes[node_ancestor]->get_number_CNV()==0) new_n_nodes_with_CNV_and_multiple_children--;

        // Hastings Ratio
        std::vector<int> nodes_containing_CNV{};
        for (int n=0;n<n_nodes;n++){
            for (auto CNV_node: nodes[n]->get_CNV_events()){
                if (CNV==CNV_node) nodes_containing_CNV.push_back(n);
            }
        }

        hastings_ratio= 1.0 / nodes_with_CNV_and_multiple_children.size() * (duplicate_CNVs.size() + nodes_with_CNV_and_multiple_children.size());
        hastings_ratio*= 
        hastings_ratio*= 1.0 * nodes_with_CNV_and_multiple_children.size() * (nodes[node_ancestor]->get_number_CNV()+1.0)
                             * children[node_ancestor].size() *descendants1.size() * descendants2.size();
        //Reverse
        hastings_ratio*= 1.0 *(duplicate_CNVs.size()+1.0) / (duplicate_CNVs.size()+1 + new_n_nodes_with_CNV_and_multiple_children);
        hastings_ratio*=  1.0 / (duplicate_CNVs.size()+1.0) / std::exp(cache_scores->log_n_choose_k(nodes_containing_CNV.size(),2)) /2.0;
    }
}

void Tree::exchange_CNV_CNLOH(){
    std::vector<int> nodes_with_event{};
    for (int i=0;i<n_nodes;i++){
        if (nodes[i]->get_number_CN_losses()>0 || nodes[i]->get_number_CNLOH()>0) nodes_with_event.push_back(i);
    }
    if (nodes_with_event.size()==0) hastings_ratio=0.0;
    else{
        int node = nodes_with_event[std::rand() % nodes_with_event.size()];
        hastings_ratio = nodes[node]->exchange_CNV_CNLOH(candidate_regions);
    }
}

void Tree::change_alleles_CNV(){
    hastings_ratio=1.0;
    std::vector<int> nodes_with_event{};
    for (int i=0;i<n_nodes;i++){
        if (nodes[i]->get_number_CNV_mut()>0) nodes_with_event.push_back(i);
    }
    if (nodes_with_event.size()==0) hastings_ratio=0.0;
    else{
        int node = nodes_with_event[std::rand() % nodes_with_event.size()];
        nodes[node]->change_alleles_CNV();
    }
}

double Tree::get_regionprobs_variance(){
    double mean = 0;
    int n_reliable_regions=0;
    for (int k=0;k<n_regions;k++){
        if (data.region_is_reliable[k]){
            n_reliable_regions++;
            mean+=region_probabilities[k] * n_regions;
        }
    }
    mean = mean / n_reliable_regions;
    std::cout<<"mean "<<mean<<std::endl;
    double variance=0;
    for (int k=0;k<n_regions;k++){
        if (data.region_is_reliable[k]){
            std::cout<<k<<": "<<region_probabilities[k] * n_regions<<std::endl;
            variance+= std::pow(region_probabilities[k]*n_regions - mean ,2);
        }
    }
    variance = variance/n_reliable_regions;
    return variance;
}

