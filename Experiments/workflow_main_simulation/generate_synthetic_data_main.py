import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import nbinom, betabinom


def generate_data(basename,method,n_nodes,n_cells,n_regions, n_SNVs,n_CNVs, n_CNLOHs, \
                    nodeprobs_concentration,dropout_rate_avg,dropout_rate_sigma,\
                    sequencing_depth,theta=10000,regionprob_sigma=0,doublet_rate=0.0,seed=0):
    np.random.seed(seed)
    sequencing_error_rate = 0.01
    # Generate the tree structure
    prufer_code = np.random.choice(n_nodes,size=n_nodes-2,replace=True)
    parents = [-1] * n_nodes
    nodes_to_add = set(range(n_nodes-1))
    for i in range(len(prufer_code)):
        n=0
        while (not n in nodes_to_add) or (n in prufer_code[i:]): n+=1
        parents[n] = prufer_code[i]
        nodes_to_add.remove(n)
    
    parents[list(nodes_to_add)[0]] = n_nodes-1
    # set the root to be 0
    parents[0],parents[n_nodes-1] = parents[n_nodes-1],parents[0]
    for i in range(len(parents)):
        if parents[i]==0: parents[i]= n_nodes-1
        elif parents[i]==n_nodes-1: parents[i]=0
    children = [[] for x in range(n_nodes)]
    for i in range(1,n_nodes): children[parents[i]].append(i)
    DFT_order=[]
    stack = [0]
    while stack!=[]:
        top = stack.pop()
        DFT_order.append(top)
        for child in children[top]: stack.append(child)


    # Assign SNVs to regions
    SNV_to_region = sorted(list(np.random.randint(0,n_regions,n_SNVs)))
    region_to_SNVs = [ [] for i in range(n_regions)]
    for i in range(n_SNVs):
        region_to_SNVs[SNV_to_region[i]].append(i)


    # Add somatic events to the nodes
    SNVs = [[] for x in range(n_nodes)]
    CNVs = [[] for x in range(n_nodes)] # each CNV is a triplet (region,gain_loss,alleles)
    CNLOHs = [[] for x in range(n_nodes)] # each CNLOH is a pair (region,alleles)

    # Assign events to nodes
    # Each node must have at least one event (SNV or CNV)
    nodes = list(range(1,n_nodes)) + list(np.random.choice(range(1,n_nodes),n_SNVs+n_CNVs-n_nodes+1,replace=True))
    np.random.shuffle(nodes)
    
    # SNVs
    for i in range(n_SNVs):
        SNVs[nodes[i]].append(i)
    for n in range(n_nodes):
        SNVs[n] = sorted(SNVs[n])

    def mut_in_ancestors(SNV,node):
        if node < 0: return False
        elif SNV in SNVs[node]: return True
        elif node<=0: return False
        else: return mut_in_ancestors(SNV,parents[node])

    # CNVs
    regions_with_CNV = sorted(np.random.choice(n_regions,n_CNVs,replace=False))
    for i in range(n_CNVs):
        n = nodes[n_SNVs+i]
        loss_gain = np.random.choice([-1,1])
        region = regions_with_CNV[i]
        alleles = []
        # CNVs can only affect alleles present in this node.
        for snv in region_to_SNVs[region]:
            if mut_in_ancestors(snv,parents[n]): alleles.append(np.random.randint(0,2))
            else: alleles.append(0)

        CNVs[n].append( (region , loss_gain , alleles) )

    # CNLOH can only affect regions containing SNVs.
    regions_with_SNV = []
    for k in range(n_regions):
        if len(region_to_SNVs[k])>0 and (not k in regions_with_CNV): regions_with_SNV.append(k)

    n_CNLOHs = min(n_CNLOHs,len(regions_with_SNV))
    
    regions_with_CNLOH = sorted(np.random.choice(regions_with_SNV,n_CNLOHs,replace=False))
    for i in range(n_CNLOHs):
        region = regions_with_CNLOH[i]
        possible_nodes = [] # Can only add the CNLOH event to a node that has a mutated allele for one of the SNVs in the region
        for n in range(1,n_nodes):
            OneSNV_in_ancestors = False
            for mut in region_to_SNVs[region]:
                OneSNV_in_ancestors = OneSNV_in_ancestors or mut_in_ancestors(mut,n)
            if OneSNV_in_ancestors: possible_nodes.append(n)
        if len(possible_nodes)>0:
            node = np.random.choice(possible_nodes)
            alleles = []
            for snv in region_to_SNVs[region]:
                 # SNV followed by a CNLOH in the next node can also be explained by two parallel branches without CNLOH, so we exclude this case
                if mut_in_ancestors(snv,parents[node]) and len(children[parents[node]])>1: alleles.append(np.random.randint(0,2))
                elif mut_in_ancestors(snv,parents[parents[node]]): alleles.append(np.random.randint(0,2))
                else: alleles.append(0)
            CNLOHs[node].append((region,alleles))

    # Compute the genotypes of each node
    cn_regions = 2*np.ones((n_regions,n_nodes),dtype=int)
    n_ref_allele = 2*np.ones((n_SNVs,n_nodes),dtype=int)
    n_alt_allele =  np.zeros((n_SNVs,n_nodes),dtype=int)
    for n in DFT_order:
        if n!=0:
            cn_regions[:,n] = cn_regions[:,parents[n]] 
            n_ref_allele[:,n] = n_ref_allele[:,parents[n]]
            n_alt_allele[:,n] = n_alt_allele[:,parents[n]]
        for mut in SNVs[n]:
            if n_ref_allele[mut,n]>0:
                n_ref_allele[mut,n]-=1
                n_alt_allele[mut,n]+=1
        for CNV in CNVs[n]:
            region,gain_loss,alleles = CNV
            cn_regions[region,n] +=gain_loss
            for i in range(len(alleles)):
                if alleles[i]==0:
                    n_ref_allele[region_to_SNVs[region][i],n]+=gain_loss
                else:
                    if n_alt_allele[region_to_SNVs[region][i],n]==0:
                        raise Exception("Invalid CNV: affects the alt allele, but its copy number before CNV is 0")
                    n_alt_allele[region_to_SNVs[region][i],n]+=gain_loss


            if len(alleles)>0 and (n_ref_allele[region_to_SNVs[region][0],n]<0 or n_alt_allele[region_to_SNVs[region][0],n]<0):
                raise Exception("Invalid CNV: copy number of one allele < 0.")
            if len(alleles)>0 and (n_ref_allele[region_to_SNVs[region][0],n]+n_alt_allele[region_to_SNVs[region][0],n]!=cn_regions[region,n]):
                raise Exception("Copy number of region does not match copy number of locus.")
        
        for CNLOH in CNLOHs[n]:
            region,alleles = CNLOH
            for i in range(len(alleles)):
                if alleles[i]==0:
                    if n_alt_allele[region_to_SNVs[region][i],n]>0: #if lose the ref allele but do not have the alt allele, do not do anything.
                        n_ref_allele[region_to_SNVs[region][i],n]-=1
                        n_alt_allele[region_to_SNVs[region][i],n]+=1
                else:
                    if n_alt_allele[region_to_SNVs[region][i],n]==0:
                        raise Exception("Invalid CNLOH: node is homozygous ref for this variant")
                    else:
                        n_ref_allele[region_to_SNVs[region][i],n]+=1
                        n_alt_allele[region_to_SNVs[region][i],n]-=1



    # Sample node probabilities
    node_probabilities = np.random.dirichlet([nodeprobs_concentration]*n_nodes) 
    # Sample dropout rates
    dropout_rates = np.abs(np.random.normal(dropout_rate_avg,dropout_rate_sigma,n_SNVs))
    # Sample region probabilities
    region_probabilities = np.abs(np.random.normal(1,regionprob_sigma,n_regions))
    region_probabilities = [max(0.3,x) for x in region_probabilities]
    # Assign cells to nodes and generate read counts
    depths = np.zeros((n_regions,n_cells),dtype=int)
    ref_reads = np.zeros((n_SNVs,n_cells),dtype=int)
    alt_reads = np.zeros((n_SNVs,n_cells),dtype=int)
    genotypes = np.zeros((n_SNVs,n_cells),dtype=int)
    for j in range(n_cells):
        if np.random.random() < doublet_rate:
            # Cell is a doublet
            node1 = int(np.random.choice(n_nodes, p=node_probabilities))
            node2 = int(np.random.choice(n_nodes, p=node_probabilities)) 
            cn_regions_cell = ( np.copy(cn_regions[:,node1]) + np.copy(cn_regions[:,node2]) ) /2
            n_ref_allele_cell = np.copy(n_ref_allele[:,node1]) + np.copy(n_ref_allele[:,node2])
            n_alt_allele_cell = np.copy(n_alt_allele[:,node1]) + np.copy(n_alt_allele[:,node2])

        else:
            # Cell is a singlet
            node = int(np.random.choice(n_nodes, p=node_probabilities))
            cn_regions_cell = np.copy(cn_regions[:,node])
            n_ref_allele_cell = np.copy(n_ref_allele[:,node]) 
            n_alt_allele_cell = np.copy(n_alt_allele[:,node])

        node = int(np.random.choice(n_nodes, p=node_probabilities))
        for k in range(n_regions):
            expected_depth = sequencing_depth * region_probabilities[k] * cn_regions_cell[k] / 2.0
            rate = np.random.gamma(shape=theta,scale=expected_depth/theta)
            depths[k,j] = np.random.poisson(rate)

        for i in range(n_SNVs):
            depth = depths[SNV_to_region[i],j]
            c_r=0 # number of copies of the ref allele that did not get dropped out
            c_a=0 # number of copies of the alt allele that did not get dropped out
            for x in range(n_ref_allele_cell[i]):
                if np.random.rand()>=dropout_rates[i]: c_r+=1
            for x in range(n_alt_allele_cell[i]):
                if (np.random.rand())>= dropout_rates[i]: c_a+=1
            if (c_r==0 and c_a==0):
                if n_ref_allele_cell[i]>0 and n_alt_allele_cell[i]:
                    if np.random.choice(2)==0: c_r+=1
                    else: c_a+=1
                elif n_ref_allele_cell[i]>0: c_r+=1
                else: c_a+=1
            f = c_a / (c_r+c_a) * (1-sequencing_error_rate) + c_r/(c_r+c_a) * sequencing_error_rate
            omega=1000 # high omega --> beta binomial ~= binomial
            a = f*omega
            b = (1-f)*omega
            alt_reads[i,j] = betabinom.rvs(depth,a,b)
            ref_reads[i,j] = depth-alt_reads[i,j]

            if ref_reads[i,j] + alt_reads[i,j]<8:
                genotypes[i,j]=3
            elif min(ref_reads[i,j],alt_reads[i,j])>0.17 * (ref_reads[i,j] + alt_reads[i,j]): 
                genotypes[i,j]=1
            elif alt_reads[i,j]>0.88 * (alt_reads[i,j] + ref_reads[i,j]):
                genotypes[i,j]=2
            elif ref_reads[i,j]>0.88 * (alt_reads[i,j] + ref_reads[i,j]):
                genotypes[i,j]=0
            else:
                genotypes[i,j]=3



    # Create files
    directory = os.path.dirname(basename)    
    if directory!="" and (not os.path.exists(directory)):
        os.makedirs(directory)

    

    # Input for COMPASS
    if method=="COMPASS":
        d={}
        d["CHR"] = [str(SNV_to_region[snv]) for snv in range(n_SNVs)]
        d["REGION"] = ["Region"+str(SNV_to_region[snv]) for snv in range(n_SNVs)]
        d["NAME"] = ["SNV"+str(snv) for snv in range(n_SNVs)]
        for j in range(n_cells):
            d[j]=[]
        for i in range(n_SNVs):
            for j in range(n_cells):
                d[j].append(str(ref_reads[i,j])+":"+str(alt_reads[i,j]))
        df_variants = pd.DataFrame(d)
        df_variants.to_csv(basename+"_variants.csv",index=False)
        
        regions = np.array(depths)
        df_regions = pd.DataFrame(regions,index = [str(i)+"_Region"+str(i) for i in range(n_regions)]).astype(int)
        df_regions.to_csv(basename+"_regions.csv",index=True,header=False)

    # Input for BiTSC2
    if method=="BITSC2":
        d={}
        d["CHR"] = [str(SNV_to_region[snv]) for snv in range(n_SNVs)]
        d["REGION"] = ["Region"+str(SNV_to_region[snv]) for snv in range(n_SNVs)]
        d["NAME"] = ["SNV"+str(snv) for snv in range(n_SNVs)]
        for j in range(n_cells):
            d[j]=[]
        for i in range(n_SNVs):
            for j in range(n_cells):
                d[j].append(str(ref_reads[i,j])+":"+str(alt_reads[i,j]))
        df_variants = pd.DataFrame(d)
        df_variants.to_csv(basename+"_variants.csv",index=False)


        DP=[]
        for i in range(n_SNVs):
            l=[]
            for j in range(n_cells):
                l.append(ref_reads[i,j]+alt_reads[i,j])
            DP.append(l)
        DP = np.array(DP)
        np.savetxt(basename+"_DP.csv",DP.astype(int),delimiter=",",fmt='%i')
        np.savetxt(basename+"_AD.csv",alt_reads.astype(int),delimiter=",",fmt='%i')

        # Create genomic segments for BITSC2
        segments = []
        start_region = 0
        end_region = 0
        while end_region < n_SNVs:
            if end_region==n_SNVs-1 or SNV_to_region[start_region] !=SNV_to_region[end_region+1]:
                segments.append((start_region+1,end_region+1))
                start_region = end_region+1
                end_region = start_region
            else:
                end_region+=1
        np.savetxt(basename+"_segments.csv",np.array(segments).astype(int),delimiter=",",fmt='%i')

    #Input for SCITE
    if method=="SCITE":
        np.savetxt(basename+"_genotypes.csv",genotypes.astype(int),delimiter=" ",fmt='%i')






    # Save tree structure, with dropout rates and node probabilities
    #parents_array = np.array(parents)
    #np.savetxt(basename+"_tree.csv",parents_array,delimiter=",",fmt='%i')
    colors = ["lightcoral","skyblue3","sandybrown","paleturquoise3","thistle","darkolivegreen3","lightpink","mediumpurple","darkseagreen3","navajowhite","gold"]
    with open(basename+"_treeTRUE.gv", "w") as file:
        file.write("digraph G{\n")
        file.write("node [color=dimgray fontsize=24 fontcolor=black fontname=Helvetica penwidth=5];\n")
        for n in range(1,n_nodes):
            file.write(str(parents[n])+" -> "+str(n)+" [color=dimgray penwidth=4 weight=2];\n")
        for n in range(n_nodes):
            label = str(n)+"[label=<"
            for mut in SNVs[n]:
                label+=str(mut)+":SNV"+str(mut)+"<br/>"
            for CNV in CNVs[n]:
                label+="CNV"
                if CNV[1]>0: label+="+1"
                else: label+="-1"
                label+= " "+str(CNV[0])+":"+  ",".join([str(x) for x in CNV[2]]) + "<br/>"
            for CNLOH in CNLOHs[n]:
                label+="CNLOH " + str(CNLOH[0])+":" + ",".join([str(x) for x in CNLOH[1]]) + "<br/>"

            if len(SNVs[n])==0 and len(CNVs[n])==0 and len(CNLOHs[n])==0: label+=" "
            label+=">];\n"
            file.write(label)
        for n in range(n_nodes):
            file.write(str(n)+" -> "+str(n+n_nodes)+ " [dir=none style=dashed weight=1 penwidth=5 color="+colors[n%len(colors)]+"];\n")
        for n in range(n_nodes):
            size = np.sqrt(100.0*node_probabilities[n]) /3.0
            file.write(str(n+n_nodes)+"[label=\"" + ("%.2f" % (100*node_probabilities[n]))+ "\\%\" style = filled width="+ str(size) \
                                    +" height="+str(size)+" color="+colors[n%len(colors)]+"];\n")
        file.write("}")


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', type = str, help='Output basename')
    parser.add_argument('--seed', type = int, default=0, help='Random seed')
    parser.add_argument('--nnodes', type = int,default = 4, help='Number of nodes in the tree')
    parser.add_argument('--ncells', type = int,default = 2000, help='Number of cells')
    parser.add_argument('--nregions', type = int,default = 50, help='Number of regions')
    parser.add_argument('--nSNVs', type = int,default = 6, help='Number of SNVs')
    parser.add_argument('--nCNVs', type = int,default = 2, help='Number of CNVs')
    parser.add_argument('--nCNLOHs', type = int,default = 0, help='Number of CNLOHs')
    parser.add_argument('--dropoutsigma', type = float,default = 0.03, help='Standard deviation when sampling the dropout rates')
    parser.add_argument('--regionprobsigma', type = float,default = 0.00, help='Standard deviation when sampling the region probabilities')
    parser.add_argument('--theta', type = float,default = 100, help='Inverse overdispersion parameter for the negative binomial (when sampling the sequencing depth)')
    parser.add_argument('--depth', type = int,default = 20, help='Sequencing depth')
    parser.add_argument('--nodeprobconcentration', type=float, default=10,help='Concentration parameter of the dirichlet distribution when sampling the node probabilities.')
    parser.add_argument('--doublet', type = float,default = 0.0, help='Doublet rate (0.0 for no doublets)')
    parser.add_argument('-m', type = str, help='Method for which the output is designed (COMPASS,BITSC2,SCITE)')
    args = parser.parse_args()
    generate_data(args.o,n_nodes=args.nnodes,n_cells=args.ncells,n_regions=args.nregions,n_SNVs=args.nSNVs,n_CNVs=args.nCNVs,n_CNLOHs=args.nCNLOHs,nodeprobs_concentration=args.nodeprobconcentration,\
        dropout_rate_avg=0.04,dropout_rate_sigma=args.dropoutsigma,seed=args.seed,sequencing_depth=args.depth,theta=args.theta,regionprob_sigma=args.regionprobsigma,doublet_rate=args.doublet,method=args.m)