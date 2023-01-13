#include <string>
#include <vector>
#include <set>
#include <random>
#include <iostream>

#include "Node.h"
#include "Structures.h"
#include "Scores.h"

//global variables
extern int n_cells;
extern int n_loci;
extern int n_regions;
extern std::vector<Cell> cells;
extern Data data;
extern Params parameters;


Node::Node(Scores* cache){
    cache_scores = cache;
    init_structures();
}

Node::Node(Node& source){
    //copy constructor
    mutations = source.mutations;
    CNA_events = source.CNA_events;
    
    n_ref_allele = source.n_ref_allele;
    n_alt_allele = source.n_alt_allele;
    cn_regions = source.cn_regions;

    affected_loci = source.affected_loci;
    affected_regions = source.affected_regions;
    attachment_scores = source.attachment_scores;
    attachment_scores_SNV = source.attachment_scores_SNV;
    attachment_scores_CNA = source.attachment_scores_CNA;

    cache_scores = source.cache_scores;

    
}

Node::Node(Node& source1, Node& source2){
    // Create a doublet from 2 nodes
    init_structures();
    cache_scores = source1.cache_scores;

    for (int i=0;i<n_loci;i++){
        // number of allele copies is the sum of the number of allele copies of the 2 nodes.
        n_ref_allele[i] = source1.n_ref_allele[i] + source2.n_ref_allele[i];
        n_alt_allele[i] = source1.n_alt_allele[i] + source2.n_alt_allele[i];

         // Actually works better if ignore the additional copy numbers, because the dropout of two copies are not independent.
       /* if (source1.n_ref_allele[i]>source2.n_ref_allele[i]) n_ref_allele[i] = source1.n_ref_allele[i];
        else n_ref_allele[i] = source2.n_ref_allele[i];
        if (source1.n_alt_allele[i]>source2.n_alt_allele[i]) n_alt_allele[i] = source1.n_alt_allele[i];
        else n_alt_allele[i] = source2.n_alt_allele[i];*/
        if (source1.n_ref_allele[i]>0 || source2.n_ref_allele[i]>0) n_ref_allele[i] = 1;
        if (source1.n_alt_allele[i]>0 || source2.n_alt_allele[i]>0) n_alt_allele[i] = 1;
        //if (n_ref_allele[i]>1) n_ref_allele[i]=1;
        //if (n_alt_allele[i]>1) n_alt_allele[i]=1;
    }
    affected_loci = source2.affected_loci;
    for (int i=0;i<n_regions;i++){
        cn_regions[i] = (source1.cn_regions[i]+source2.cn_regions[i])/2;
    }
    affected_regions = source2.affected_regions;
}


void Node::init_structures(){
    n_ref_allele.resize(n_loci);
    n_alt_allele.resize(n_loci);
    cn_regions.resize(n_regions);

    attachment_scores.resize(n_cells);
    attachment_scores_SNV.resize(n_cells);
    attachment_scores_CNA.resize(n_cells);
}

Node::~Node(){
}

void Node::update_genotype(Node* parent){

    // 1. Initialize with the parent genotype.
    if (parent==nullptr){
        // Start from diploid homozygous reference (except X chromosome for males)
        for (int i=0;i<n_loci;i++){
            if (data.locus_to_chromosome.size()>0 && data.locus_to_chromosome[i]=="X" && data.sex == "male") n_ref_allele[i]= 1;
            else n_ref_allele[i] = 2;
            n_alt_allele[i] = 0;
            
        }
        for (int k=0;k<n_regions;k++){
            if (data.region_to_chromosome.size()>0 && data.region_to_chromosome[k]=="X" && data.sex=="male") cn_regions[k]=1;
            else cn_regions[k] = 2;
        }
        
    }
    else{
        for (int i=0;i<n_loci;i++){
            n_ref_allele[i] = parent->n_ref_allele[i];
            n_alt_allele[i] = parent->n_alt_allele[i];
        }
        for (int i=0;i<n_regions;i++){
            cn_regions[i] = parent->cn_regions[i];
        }
        
    }
    // 2. Go through the events, and update the genotype accordingly.
    affected_loci.clear(); //set of loci which are different from the parent
    affected_regions.clear();

    // 2.1. Mutations
    for (const int& mutated_locus: mutations){
        //somatic mutation: go from ref to alt
        if (n_ref_allele[mutated_locus]>=1){
            n_ref_allele[mutated_locus]-=1;
            n_alt_allele[mutated_locus]+=1;
            affected_loci.insert(mutated_locus);
        }
    }


    // 2.2. CNA

    // If a CNA affects an allele not present, change the affected allele
    std::vector<std::tuple<int,int,std::vector<int>>> CNAs_to_remove{};
    std::vector<std::tuple<int,int,std::vector<int>>> CNAs_to_add{};
    for (auto CNA: CNA_events){
        int region = std::get<0>(CNA);
        bool valid=true;
        std::vector<int> alleles = std::get<2>(CNA);
        std::vector<int> new_alleles{};
        for (int i=0; i< data.region_to_loci[region].size();i++){
            if (alleles[i]==0 && n_ref_allele[data.region_to_loci[region][i]]==0){
                valid = false;
                new_alleles.push_back(1);
            }
            else if (alleles[i]==1 &&  n_alt_allele[data.region_to_loci[region][i]]==0){
                valid = false;
                new_alleles.push_back(0);
            }
            else new_alleles.push_back(alleles[i]);
        }
        if (!valid){
            CNAs_to_remove.push_back(CNA);
            CNAs_to_add.push_back(std::make_tuple(region,std::get<1>(CNA),new_alleles));
        }
    }

    for (auto CNA: CNAs_to_remove) remove_CNA(CNA);
    for (auto CNA: CNAs_to_add) add_CNA(CNA);




    all_CNA_events_valid = true;
    for (const std::tuple<int,int,std::vector<int>> CNA: CNA_events){
        int region = std::get<0>(CNA);
        int gain_loss = std::get<1>(CNA);
        const std::vector<int>& alleles = std::get<2>(CNA);
        if (parameters.verbose) std::cout<<"CNA"<<gain_loss<<" in " <<data.region_to_name[region]<<std::endl;
        if (gain_loss!=0) affected_regions.insert(region);

        // Check that the CNA is valid (region ends up with a copy number in {1,2,3} and affected alleles have copy number >0)
        bool valid_CNA=(cn_regions[region]+gain_loss>=1 & cn_regions[region]+gain_loss<=3);
        for (int i=0;i<data.region_to_loci[region].size();i++){
            int locus = data.region_to_loci[region][i];
            if (alleles[i]==0 && n_ref_allele[locus]==0) valid_CNA=false;
            if (alleles[i]==1 && n_alt_allele[locus]==0) valid_CNA=false;
        }

        // Apply the CNA
        if (valid_CNA) {
            cn_regions[region]+=gain_loss;
            for (int i=0;i<data.region_to_loci[region].size();i++){
                affected_loci.insert(data.region_to_loci[region][i]);
                int locus = data.region_to_loci[region][i];
                if (gain_loss!=0){
                    if (alleles[i]==0) n_ref_allele[locus]+=gain_loss;
                    else if (alleles[i]==1) n_alt_allele[locus]+=gain_loss;
                }
                else if (n_ref_allele[locus]>0 &&  n_alt_allele[locus]>0){ // copy neutral
                    if (alleles[i]==0) {
                        n_ref_allele[locus]-=1;
                        n_alt_allele[locus]+=1;
                    }
                    else{
                        n_ref_allele[locus]+=1;
                        n_alt_allele[locus]-=1;
                    }
                }
                
            }   
        }
        else{
            all_CNA_events_valid = false;
        }
    }
}

void Node::compute_attachment_scores(bool use_CNA,const std::vector<double>& dropout_rates_ref,
                        const std::vector<double>& dropout_rates_alt, const std::vector<double>& region_probabilities){
    // Compute the attachment score of a cell to a node, starting from scratch (for the root).
    for (int j=0;j<n_cells;j++){
        attachment_scores_SNV[j]=0.0;
        attachment_scores_CNA[j]=0.0;
    }
    
    std::vector<double> temp_scores{};
    for (int i=0; i<n_loci;i++){
        temp_scores = cache_scores->compute_SNV_loglikelihoods(n_ref_allele[i],n_alt_allele[i],i,dropout_rates_ref[i],dropout_rates_alt[i]);
        for (int j=0;j<n_cells;j++){
            attachment_scores_SNV[j]+=temp_scores[j];
        }
    }
    for (int j=0;j<n_cells;j++){
        attachment_scores[j] = attachment_scores_SNV[j];
    }

    // No CNA at the root allowed, so no need to compute CNA scores at the root (constant offset)
}

void Node::compute_attachment_scores_parent(bool use_CNA, Node* parent,const std::vector<double>& dropout_rates_ref,
                            const std::vector<double>& dropout_rates_alt, const std::vector<double>& region_probabilities,bool recompute_CNA_scores){
    // Compute the attachment score of a cell to a node, starting from the attachment score of the cell to the parent of the current node.
    // Only update the score for loci and regions where the genotype differs from the parent.
    // The CNA scores need to be computed once per tree, because they do not depend on the parameters inferred during the MCMC.
    attachment_scores_SNV = parent->attachment_scores_SNV;
    std::vector<double> temp_scores{};
    for (int i: affected_loci){
        temp_scores = cache_scores->compute_SNV_loglikelihoods(parent->n_ref_allele[i],parent->n_alt_allele[i],i,dropout_rates_ref[i],dropout_rates_alt[i]);
        for (int j=0;j<n_cells;j++) attachment_scores_SNV[j]-=temp_scores[j];
        temp_scores = cache_scores->compute_SNV_loglikelihoods(n_ref_allele[i],n_alt_allele[i],i,dropout_rates_ref[i],dropout_rates_alt[i]);                                                
        for (int j=0;j<n_cells;j++) attachment_scores_SNV[j]+=temp_scores[j];
    }
   
    if (use_CNA){
        if (true){
            if (recompute_CNA_scores){ // Only compute the CNA scores once per tree.
                attachment_scores_CNA = parent->attachment_scores_CNA;
                if (affected_regions.size()>0){ // If there are no CNA, can keep the CNA score of the parent.
                    double normalization_factor=0.0;
                    double normalization_factor_parent=0.0;
                    for (int k=0;k<n_regions;k++){
                        if (data.region_is_reliable[k]){
                            normalization_factor+= region_probabilities[k] *cn_regions[k];
                            normalization_factor_parent+=region_probabilities[k] * parent->cn_regions[k];
                        }
                    }
                    for (int k=0;k<n_regions;k++){
                        if (data.region_is_reliable[k]){
                            temp_scores = cache_scores->compute_CNA_loglikelihoods(k,region_probabilities[k] * parent->cn_regions[k]/ normalization_factor_parent);
                            for (int j=0;j<n_cells;j++) attachment_scores_CNA[j]-=temp_scores[j];
                            temp_scores = cache_scores->compute_CNA_loglikelihoods(k,region_probabilities[k] * cn_regions[k]/normalization_factor);
                            for (int j=0;j<n_cells;j++) attachment_scores_CNA[j]+=temp_scores[j];
                        }
                    }
                }
            }
        }
        else{
            for (int k: affected_regions){
                if (data.region_is_reliable[k]){
                    temp_scores = cache_scores->compute_CNA_loglikelihoods(k,region_probabilities[k] * parent->cn_regions[k]/2.0);
                    for (int j=0;j<n_cells;j++) attachment_scores_CNA[j]-=temp_scores[j];
                    temp_scores = cache_scores->compute_CNA_loglikelihoods(k,region_probabilities[k] * cn_regions[k]/2.0);
                    for (int j=0;j<n_cells;j++) attachment_scores_CNA[j]+=temp_scores[j];
                }
            }
        }
        for (int j=0;j<n_cells;j++){
            attachment_scores[j] = attachment_scores_SNV[j] + attachment_scores_CNA[j];
        }
    }
    else{
        for (int j=0;j<n_cells;j++){
            attachment_scores[j] = attachment_scores_SNV[j];
        }
    }
    
}


int Node::remove_random_mutation(){
    // Randomly remove one of the mutations and return it
    // This method should only be called if the node has at least one somatic mutation.
    int idx = std::rand()%mutations.size();
    int mutation = mutations[idx];
    mutations.erase(mutations.begin()+idx);
    return mutation;
}


std::tuple<int,int,std::vector<int>> Node::remove_random_CNA(){
    // Randomly remove one of the existing CNA events and return it. 
    // This method should only be called if the node has at least one CNA event.
    int index_to_remove = std::rand()%CNA_events.size();
    std::tuple<int,int,std::vector<int>> CNA = CNA_events[index_to_remove];
    CNA_events.erase(CNA_events.begin()+index_to_remove); // remove only one occurence of the event (multiset)
    return CNA;
}


double Node::exchange_Loss_CNLOH(std::vector<int> candidate_regions){
    double hastings_ratio=1.0;
    
    std::multiset<std::tuple<int,int,std::vector<int>>> CNLOH_events_exchangeable{}; // CNLOH events in a region that is allowed to contain a CNA
    for (auto CNA: CNA_events){
        if (std::get<1>(CNA)==0){
            for (int k: candidate_regions){
                if (std::get<0>(CNA)==k) CNLOH_events_exchangeable.insert(CNA);
            }
        }
    }
    int n_CNLOH = CNLOH_events_exchangeable.size();

    std::multiset<std::tuple<int,int,std::vector<int>>> CN_losses{}; // can only transform a copy number loss (with some loci) into a CNLOH
    for (auto CNA: CNA_events){
        if (std::get<1>(CNA)==-1 && std::get<2>(CNA).size()>0) CN_losses.insert(CNA);
    }
    int n_CN_loss = CN_losses.size();
    if (n_CN_loss + n_CNLOH ==0){
        return 0.0;
    }
    int event_ind = std::rand()%(n_CN_loss+n_CNLOH);
    if (event_ind < n_CNLOH){ // transform CNLOH event into a CNA event
        std::tuple<int,int,std::vector<int>> CNLOH_event = *std::next(CNLOH_events_exchangeable.begin(), event_ind);
        int region = std::get<0>(CNLOH_event);
        std::vector<int> lost_alleles = std::get<2>(CNLOH_event);

        // Replace the CNLOH with a loss
        for (int i=0;i<CNA_events.size();i++){
            if (std::get<0>(CNA_events[i])==region){
                CNA_events[i] = std::make_tuple(region,-1,lost_alleles);
            }
        }
    }
    else{ // transform Loss into a CNLOH event
        std::tuple<int,int,std::vector<int>> CNA_event = *std::next(CN_losses.begin(), event_ind-n_CNLOH);
        int region = std::get<0>(CNA_event);
        std::vector<int> alleles = std::get<2>(CNA_event);
        // Replace the loss with a CNLOH
        for (int i=0;i<CNA_events.size();i++){
            if (std::get<0>(CNA_events[i])==region){
                CNA_events[i] = std::make_tuple(region,0,alleles);
            }
        }
    }
    return hastings_ratio;
}

void Node::change_alleles_CNA(){
    // For a CNA event, change which alleles are lost or gained.
    
    // Choose one CNA event affecting a region which contains variants
    std::vector<std::tuple<int,int,std::vector<int>>> CNA_with_muts{};
    for (auto CNA:CNA_events){
        if (std::get<2>(CNA).size()>0) CNA_with_muts.push_back(CNA);
    }
    std::tuple<int,int,std::vector<int>> CNA = CNA_with_muts[std::rand()%CNA_with_muts.size()];
    // Replace the affected alleles
    int region = std::get<0>(CNA);
    int gain_loss = std::get<1>(CNA);
    std::vector<int> alleles{};
    for (int i=0;i<data.region_to_loci[region].size();i++){
        int allele = std::rand()%2;
        alleles.push_back(allele);
    }
    for (int i=0;i<CNA_events.size();i++){
        if (std::get<0>(CNA_events[i])==region){
            CNA_events[i] = std::make_tuple(region,gain_loss,alleles);
        }
    }
}

void Node::change_alleles_CNA_locus(int locus, bool heterozygous){
    // Change the allele affected by a CNA, only at a particular locus
    int region = data.locus_to_region[locus];
    for (int i=0;i<CNA_events.size();i++){
        if (std::get<0>(CNA_events[i])==region){
            std::vector<int> alleles = std::get<2>(CNA_events[i]);
            std::vector<int> new_alleles{};
            for (int a=0;a<alleles.size();a++){
                if (locus==data.region_to_loci[region][a]){
                    if (heterozygous){
                        new_alleles.push_back(std::rand()%2);
                    }
                    else{
                        new_alleles.push_back(0);
                    }
                }
                else new_alleles.push_back(alleles[a]);
            }
            CNA_events[i] = std::make_tuple(region,std::get<1>(CNA_events[i]),new_alleles);
        }
    }
}




int Node::get_number_disjoint_CNA(std::vector<int> regions_successor){
    // Several events count as one if:
    // - they are adjacent (in the list of regions)
    // - they are on the same chromosome
    // - they are of the same type (Gain, Loss or CNLOH)
    int count=0;
    int last_region=-10;
    int last_type=-10;
    for (std::tuple<int,int,std::vector<int>> CNA:CNA_events){
        int region = std::get<0>(CNA);
        int type = std::get<1>(CNA);
        bool regions_adjacent = data.region_to_chromosome[region]==data.region_to_chromosome[last_region];
        if (regions_adjacent){
            if (type==0){
                regions_adjacent = (region==last_region+1);
            }
            else{
                regions_adjacent = (region==regions_successor[last_region]);
            }
        }
        if ( last_type!=-10 && ((!regions_adjacent) || type!=last_type) ){
            count++;
        }
        last_region=region;
        last_type = type;
    }
    return count;
}

int Node::get_number_disjoint_LOH(std::vector<int> regions_successor){
    // Several events count as one if:
    // - they are adjacent (in the list of regions)
    // - they are on the same chromosome
    // - they are of the same type (Gain, Loss or CNLOH)
    // Here, only count Losses and CNLOH in regions which contain a variant (resulting in a LOH)
    int count=0;
    bool segment_contains_LOH = false;
    int last_region=-10;
    int last_type=-10;
    for (std::tuple<int,int,std::vector<int>> CNA:CNA_events){
        int region = std::get<0>(CNA);
        int type = std::get<1>(CNA);
        bool regions_adjacent = data.region_to_chromosome[region]==data.region_to_chromosome[last_region];
        if (regions_adjacent){
            if (type==0){
                regions_adjacent = (region==last_region+1);
            }
            else{
                regions_adjacent = (region==regions_successor[last_region]);
            }
        }
        if ( last_type!=-10 && ( (!regions_adjacent) || type!=last_type) ){
            if (segment_contains_LOH) count++;
           segment_contains_LOH = false;
        }
        segment_contains_LOH = segment_contains_LOH || (data.region_to_loci[region].size()>0 && type<=0);
        last_region=region;
        last_type = type;
    }
    if (segment_contains_LOH) count++;
    return count;
}




std::string Node::get_label(){
    // This label (meant for graphviz) contains the list of mutations in the node.
    std::string label{};
    for (int i=0;i<n_loci;i++){
        for (int mut: mutations){
            if (i==mut){
                //check if this mutation is important or not
                bool nonsyn_mut=false;
                for (int pos=0;pos<data.locus_to_name[i].size()-1;pos++){
                    if (data.locus_to_name[i].substr(pos,2)=="p.") nonsyn_mut=true;
                }
                if (data.locus_to_name[i]=="FLT3-ITD") nonsyn_mut=true;
                if (data.locus_to_name[i].size()>12 && data.locus_to_name[i].substr(data.locus_to_name[i].size()-12,12)=="splice-donor") nonsyn_mut = true;
                bool somatic_nonsyn_mut = nonsyn_mut && data.locus_to_freq[i]==0.0;
                if (somatic_nonsyn_mut) label+="<B>";
                label+= std::to_string(mut) + ": " + data.locus_to_name[i] +"(chr"+data.locus_to_chromosome[mut]+")";
                if (somatic_nonsyn_mut) label+="</B>";
                label+="<br/>";
            }
        }
    }
    for (std::tuple<int,int,std::vector<int>> CNA:CNA_events){
        int region = std::get<0>(CNA);
        std::string type{};
        if (std::get<1>(CNA)==1) type="Gain";
        else if (std::get<1>(CNA)==-1) type="Loss";
        else type="CNLOH";
        label+= "<B>" + type + " "+std::to_string(region)+":"+ data.region_to_name[region]+ "(chr"+data.region_to_chromosome[region]+ "):";
        std::vector<int> alleles = std::get<2>(CNA);
        for (int i = 0; i<alleles.size();i++){
            label+=std::to_string(alleles[i]);
            if (i+1<alleles.size()) label+=",";
        }
        label+="</B><br/>";
    }
    if (label=="") label = " ";
    return label;
}

std::string Node::get_label_simple(std::set<int> excluded_mutations){
    // This label (meant for graphviz) contains the list of mutations in the node.
    std::string label{};
    for (int i=0;i<n_loci;i++){
        if (excluded_mutations.count(i)) continue;
        for (int mut: mutations){
            if (i==mut){
                //check if this mutation is important or not
                bool nonsyn_mut=false;
                for (int pos=0;pos<data.locus_to_name[i].size()-1;pos++){
                    if (data.locus_to_name[i].substr(pos,2)=="p.") nonsyn_mut=true;
                }
                if (data.locus_to_name[i]=="FLT3-ITD") nonsyn_mut=true;
                if (data.locus_to_name[i].size()>12 && data.locus_to_name[i].substr(data.locus_to_name[i].size()-12,12)=="splice-donor") nonsyn_mut = true;
                bool somatic_nonsyn_mut = nonsyn_mut && data.locus_to_freq[i]==0.0;
                if (somatic_nonsyn_mut) label+="<B>";
                label+= data.locus_to_name[i]+"(chr"+data.locus_to_chromosome[mut]+")";
                if (somatic_nonsyn_mut) label+="</B>";
                label+="<br/>";
            }
        }
    }
    for (std::tuple<int,int,std::vector<int>> CNA:CNA_events){
        int region = std::get<0>(CNA);
        std::string type{};
        if (std::get<1>(CNA)==1) type="Gain";
        else if (std::get<1>(CNA)==-1) type="Loss";
        else type="CNLOH";
        label+= "<B>" + type + " "+ data.region_to_name[region];
        std::vector<int> alleles = std::get<2>(CNA);
        if (alleles.size()>0) label+=":";
        for (int i = 0; i<alleles.size();i++){
            if (alleles[i]==0) label+="REF";
            else label+="ALT";
            if (i+1<alleles.size()) label+=",";
        }
        label+=" (chr" + data.region_to_chromosome[region]+")";
        label+="</B><br/>";
    }
    if (label=="") label = " ";
    return label;
}
