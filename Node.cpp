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
    CNLOH_events = source.CNLOH_events;
    CNV_events = source.CNV_events;
    
    n_ref_allele = source.n_ref_allele;
    n_alt_allele = source.n_alt_allele;
    cn_regions = source.cn_regions;

    affected_loci = source.affected_loci;
    affected_regions = source.affected_regions;
    attachment_scores = source.attachment_scores;

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
    // 2.2. CN-LOH
    for (const std::pair<int,std::vector<int>>& CNLOH: CNLOH_events){
        int region = CNLOH.first;
        for (int i=0;i<CNLOH.second.size();i++){
            int locus = data.region_to_loci[region][i];
            if (n_ref_allele[locus]>0 && n_alt_allele[locus]>0){ // must be heterozygous to have a CN-LOH event
            affected_loci.insert(locus);
            if (CNLOH.second[i]==0){ // Lose the ref allele
                n_ref_allele[locus]--;
                n_alt_allele[locus]++;
            }
            else{ // Lose the alt allele
                n_ref_allele[locus]++;
                n_alt_allele[locus]--;
            }
        }
        // don't do anything special if we have a CN-LOH event but were already homozygous.
        }
        
    }
    // 2.3. CNV

    // If a CNV affects an allele not present, change the affected allele
    std::vector<std::tuple<int,int,std::vector<int>>> CNVs_to_remove{};
    std::vector<std::tuple<int,int,std::vector<int>>> CNVs_to_add{};
    for (auto CNV: CNV_events){
        int region = std::get<0>(CNV);
        bool valid=true;
        std::vector<int> alleles = std::get<2>(CNV);
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
            CNVs_to_remove.push_back(CNV);
            CNVs_to_add.push_back(std::make_tuple(region,std::get<1>(CNV),new_alleles));
        }
    }

    for (auto CNV: CNVs_to_remove) remove_CNV(CNV);
    for (auto CNV: CNVs_to_add) add_CNV(CNV);




    all_CNV_events_valid = true;
    for (const std::tuple<int,int,std::vector<int>> CNV: CNV_events){
        int region = std::get<0>(CNV);
        int gain_loss = std::get<1>(CNV);
        const std::vector<int>& alleles = std::get<2>(CNV);
        if (parameters.verbose) std::cout<<"CNV"<<gain_loss<<" in " <<data.region_to_name[region]<<std::endl;
        affected_regions.insert(region);

        // Check that the CNV is valid (region ends up with a copy number in {1,2,3} and affected alleles have copy number >0)
        bool valid_CNV=(cn_regions[region]+gain_loss>=1 & cn_regions[region]+gain_loss<=3);
        for (int i=0;i<data.region_to_loci[region].size();i++){
            int locus = data.region_to_loci[region][i];
            if (alleles[i]==0 && n_ref_allele[locus]==0) valid_CNV=false;
            if (alleles[i]==1 && n_alt_allele[locus]==0) valid_CNV=false;
        }
        // check that the CNV is among the candidate CNVs.
        /*bool CNV_in_candidates=false;
        for (int k: candidate_regions) CNV_in_candidates = CNV_in_candidates || (region==k);
        valid_CNV = valid_CNV && CNV_in_candidates;*/

        // Apply the CNV
        if (valid_CNV) {
            cn_regions[region]+=gain_loss;
            for (int i=0;i<data.region_to_loci[region].size();i++){
                affected_loci.insert(data.region_to_loci[region][i]);
                int locus = data.region_to_loci[region][i];
                if (alleles[i]==0) n_ref_allele[locus]+=gain_loss;
                else if (alleles[i]==1) n_alt_allele[locus]+=gain_loss;
            }   
        }
        else{
            all_CNV_events_valid = false;
        }
    }
}

void Node::compute_attachment_scores(bool use_CNV,const std::vector<double>& dropout_rates_ref,
                        const std::vector<double>& dropout_rates_alt, const std::vector<double>& region_probabilities){
    // Compute the attachment score of a cell to a node, starting from scratch (for the root).
    for (int j=0;j<n_cells;j++) attachment_scores[j]=0;
    
    std::vector<double> temp_scores{};
    for (int i=0; i<n_loci;i++){
        temp_scores = cache_scores->compute_SNV_loglikelihoods(n_ref_allele[i],n_alt_allele[i],i,dropout_rates_ref[i],dropout_rates_alt[i]);
        for (int j=0;j<n_cells;j++) attachment_scores[j]+=temp_scores[j];
    }

    // No CNV at the root allowed, so now need to compute CNV scores at the root (constant offset)
}

void Node::compute_attachment_scores_parent(bool use_CNV, Node* parent,const std::vector<double>& dropout_rates_ref,
                            const std::vector<double>& dropout_rates_alt, const std::vector<double>& region_probabilities){
    // Compute the attachment score of a cell to a node, starting from the attachment score of the cell to the parent of the current node.
    // Only update the score for loci and regions where the genotype differs from the parent.
    attachment_scores = parent->attachment_scores;
    std::vector<double> temp_scores{};
    for (int i: affected_loci){
        temp_scores = cache_scores->compute_SNV_loglikelihoods(parent->n_ref_allele[i],parent->n_alt_allele[i],i,dropout_rates_ref[i],dropout_rates_alt[i]);
        for (int j=0;j<n_cells;j++) attachment_scores[j]-=temp_scores[j];
        temp_scores = cache_scores->compute_SNV_loglikelihoods(n_ref_allele[i],n_alt_allele[i],i,dropout_rates_ref[i],dropout_rates_alt[i]);
                                                        
        for (int j=0;j<n_cells;j++) attachment_scores[j]+=temp_scores[j];
    }
   
    if (use_CNV){
        for (int k: affected_regions){
            if (data.region_is_reliable[k]){
                temp_scores = cache_scores->compute_CNV_loglikelihoods(k,region_probabilities[k] * parent->cn_regions[k]/2.0);
                for (int j=0;j<n_cells;j++) attachment_scores[j]-=temp_scores[j];
                temp_scores = cache_scores->compute_CNV_loglikelihoods(k,region_probabilities[k] * cn_regions[k]/2.0);
                for (int j=0;j<n_cells;j++) attachment_scores[j]+=temp_scores[j];
            }
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

std::pair<int,std::vector<int>> Node::remove_random_CNLOH(){
    // Randomly remove one of the existing CNLOH events and return it. 
    // This method should only be called if the node has at least one CNLOH event.
    int index_to_remove = std::rand()%CNLOH_events.size();
    std::pair<int,std::vector<int>> event = *std::next(CNLOH_events.begin(), index_to_remove);
    CNLOH_events.erase(CNLOH_events.lower_bound(event)); // remove only one occurence of the event (multiset)
    return event;
}


std::tuple<int,int,std::vector<int>> Node::remove_random_CNV(){
    // Randomly remove one of the existing CNV events and return it. 
    // This method should only be called if the node has at least one CNV event.
    int index_to_remove = std::rand()%CNV_events.size();
    std::tuple<int,int,std::vector<int>> event = *std::next(CNV_events.begin(), index_to_remove);
    CNV_events.erase(CNV_events.lower_bound(event)); // remove only one occurence of the event (multiset)
    return event;
}

void Node::remove_CNVs_in_region(int region){
    // Remove the CNVs affecting this region, if there are any.
    std::vector<std::tuple<int,int,std::vector<int>>> CNVs_to_delete{};
    for (auto CNV: CNV_events){
        if (std::get<0>(CNV)==region) CNVs_to_delete.push_back(CNV);
    }
    for (auto CNV: CNVs_to_delete){
        CNV_events.erase(CNV_events.lower_bound(CNV));
    }
}

double Node::exchange_CNV_CNLOH(std::vector<int> candidate_regions){
    double hastings_ratio=1.0;
    
    std::multiset<std::pair<int,std::vector<int>>> CNLOH_events_exchangeable{}; // CNLOH events in a region that is allowed to contain a CNV
    for (auto CNLOH: CNLOH_events){
        for (int k: candidate_regions){
            if (CNLOH.first==k) CNLOH_events_exchangeable.insert(CNLOH);
        }
    }
    int n_CNLOH = CNLOH_events_exchangeable.size();

    std::multiset<std::tuple<int,int,std::vector<int>>> CN_losses{}; // can only transform a copy number loss (with some loci) into a CNLOH
    for (auto CNV: CNV_events){
        if (std::get<1>(CNV)==-1 && std::get<2>(CNV).size()>0) CN_losses.insert(CNV);
    }
    int n_CN_loss = CN_losses.size();
    if (n_CN_loss + n_CNLOH ==0){
        return 0.0;
    }
    int event_ind = std::rand()%(n_CN_loss+n_CNLOH);
    if (event_ind < n_CNLOH){ // transform CNLOH event into a CNV event
        std::pair<int,std::vector<int>> CNLOH_event = *std::next(CNLOH_events_exchangeable.begin(), event_ind);
        int region = CNLOH_event.first;
        std::vector<int> lost_alleles = CNLOH_event.second;
        // Check if there is already a CNV affecting this region. In this case, might merge both CNLOH and CNV into a different CNV.
        // This can happen when we have a CNLOH of the ref allele and a copy number loss of one alt allele instead of a copy number loss of the ref allele.
        for (auto CNV: CN_losses){
            if (std::get<0>(CNV) == region){
                hastings_ratio*=2.0;
                if (std::rand()%2==0) CNV_events.erase(CNV_events.lower_bound(CNV));
            }
        }
        CNLOH_events.erase(CNLOH_events.lower_bound(CNLOH_event));
        CNV_events.insert(std::make_tuple(region,-1,lost_alleles));
    }
    else{ // transform CNV event into a CNLOH event
        std::tuple<int,int,std::vector<int>> CNV_event = *std::next(CN_losses.begin(), event_ind-n_CNLOH);
        int region = std::get<0>(CNV_event);
        std::vector<int> alleles = std::get<2>(CNV_event);
        CNV_events.erase(CNV_events.lower_bound(CNV_event));
        CNLOH_events.insert(std::make_pair(region,alleles));
    }
    return hastings_ratio;
}

void Node::change_alleles_CNV(){
    // For a CNV event, change which alleles are lost or amplified.
    
    // Choose one CNV event affecting a region which contains variants
    std::vector<std::tuple<int,int,std::vector<int>>> CNV_with_muts{};
    for (auto CNV:CNV_events){
        if (std::get<2>(CNV).size()>0) CNV_with_muts.push_back(CNV);
    }
    std::tuple<int,int,std::vector<int>> CNV = CNV_with_muts[std::rand()%CNV_with_muts.size()];
    // Replace the affected alleles
    int region = std::get<0>(CNV);
    int gain_loss = std::get<1>(CNV);
    CNV_events.erase(CNV_events.lower_bound(CNV));
    std::vector<int> alleles{};
    for (int i=0;i<data.region_to_loci[region].size();i++){
        int allele = std::rand()%2;
        alleles.push_back(allele);
    }
    CNV_events.insert(std::make_tuple(region,gain_loss,alleles));

    // If a CNLOH affects the same region, can remove it 
    // (in case the CNV amplifies the ref allele and there is a CNLOH from ref to alt, instead of only having a CN gain of the alt allele)
    std::vector<std::pair<int,std::vector<int>>> CNLOH_to_remove{};
    for (auto CNLOH: CNLOH_events){
        if (CNLOH.first==region){
            if (std::rand()%2 ==0) CNLOH_to_remove.push_back(CNLOH);
        }
    }
    for (auto CNLOH: CNLOH_to_remove){
        CNLOH_events.erase(CNLOH_events.lower_bound(CNLOH));
    }
}



int Node::get_number_disjoint_CNLOH(){
    int count=0;
    int lastpos=-10;
    for (std::pair<int,std::vector<int>> CNLOH:CNLOH_events){
        if (data.region_to_chromosome.size()==0 || CNLOH.first!=lastpos+1 || data.region_to_chromosome[CNLOH.first]!=data.region_to_chromosome[lastpos]){
            if (data.region_to_chromosome.size()>0 && lastpos>=0 && data.region_to_chromosome[CNLOH.first]==data.region_to_chromosome[lastpos]){
                // If there are two CNLOH on the same chromosome, but not on adjacent regions, they might still be considered as one event if there are no 
                // regions with variants between them.
                bool no_region_with_variant_between=lastpos<CNLOH.first;
                for (int k=lastpos+1;k<CNLOH.first;k++){
                    if (data.region_to_loci[k].size()>0) no_region_with_variant_between=false;
                }
                if (no_region_with_variant_between) continue;
            }
            count++;
        }
        lastpos=CNLOH.first;
    }
    return count;
}

int Node::get_number_disjoint_CNV(std::vector<int> regions_successor){
    // Several events count as one if:
    // - they are adjacent (in the list of regions)
    // - they are on the same chromosome
    // - they have the same sign (both gains or both losses)
    int count=0;
    int last_region=regions_successor.size()-1;
    int last_gain_loss=-10;
    for (std::tuple<int,int,std::vector<int>> CNV:CNV_events){
        int region = std::get<0>(CNV);
        int gain_loss = std::get<1>(CNV);
        if ( (last_region<regions_successor.size() && region!=regions_successor[last_region]) || gain_loss!=last_gain_loss){
            count++;
        }
        last_region=region;
        last_gain_loss = gain_loss;
    }
    return count;
}

int Node::get_number_disjoint_CNV_LOH(std::vector<int> regions_successor){
    // Several events count as one if:
    // - they are adjacent (in the list of regions)
    // - they are on the same chromosome
    // - they have the same sign (both gains or both losses)
    // Here, only count CN losses which contain a variant (resulting in a LOH)
    int count=0;
    bool segment_contains_LOH = false;
    int last_region=regions_successor.size()-1;
    int last_gain_loss=-10;
    for (std::tuple<int,int,std::vector<int>> CNV:CNV_events){
        int region = std::get<0>(CNV);
        int gain_loss = std::get<1>(CNV);
        if ( (last_region<regions_successor.size() && region!=regions_successor[last_region]) || gain_loss!=last_gain_loss){
            if (segment_contains_LOH) count++;
           segment_contains_LOH = false;
        }
        segment_contains_LOH = segment_contains_LOH || (data.region_to_loci[region].size()>0 && gain_loss==-1);
        last_region=region;
        last_gain_loss = gain_loss;
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
    for (std::pair<int,std::vector<int>> CNLOH: CNLOH_events){
        int region = CNLOH.first;
        label+= "<B>CNLOH " + std::to_string(region)+":" +data.region_to_name[region] + "(chr"+data.region_to_chromosome[region]+ "):";
        for (int i=0;i<CNLOH.second.size();i++){
            label+=std::to_string(CNLOH.second[i]);
            if (i+1<CNLOH.second.size()) label+=",";
        }
        label += "</B><br/>";
    }
    for (std::tuple<int,int,std::vector<int>> CNV:CNV_events){
        int region = std::get<0>(CNV);
        std::string sign{};
        if (std::get<1>(CNV)==1) sign="+1";
        else sign="-1";
        label+= "<B>CNV" + sign + " "+std::to_string(region)+":"+ data.region_to_name[region]+ "(chr"+data.region_to_chromosome[region]+ "):";
        std::vector<int> alleles = std::get<2>(CNV);
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
    for (std::pair<int,std::vector<int>> CNLOH: CNLOH_events){
        int region = CNLOH.first;
        label+= "<B>CNLOH "  +data.region_to_name[region] + "(chr"+data.region_to_chromosome[region]+ "):";
        for (int i=0;i<CNLOH.second.size();i++){
            if (CNLOH.second[i]==0) label+="REF";
            else label+="ALT";
            if (i+1<CNLOH.second.size()) label+=",";
        }
        label += "</B><br/>";
    }
    for (std::tuple<int,int,std::vector<int>> CNV:CNV_events){
        int region = std::get<0>(CNV);
        std::string sign{};
        if (std::get<1>(CNV)==1) sign="+1";
        else sign="-1";
        label+= "<B>CNV" + sign + " "+ data.region_to_name[region];
        std::vector<int> alleles = std::get<2>(CNV);
        if (alleles.size()>0) label+=":";
        std::cout<<data.region_to_name[region]<<": "<<alleles.size();
        for (int i = 0; i<alleles.size();i++){
            std::cout<<alleles[i]<<std::endl;
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
