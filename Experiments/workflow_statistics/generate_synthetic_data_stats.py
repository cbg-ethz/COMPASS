import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import nbinom, betabinom


def generate_data(basename,n_nodes,n_cells,n_SNVs,n_CNVs, \
                    nodeprobs_concentration,dropout_rate_avg,dropout_rate_sigma,\
                    sequencing_depth,theta=10000,regionprob_sigma=0,seed=0):
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


    # Add somatic events to the nodes
    muts = [[] for x in range(n_nodes)]
    CNVs = [[] for x in range(n_nodes)] # each CNV is a triplet (locus,gain_loss,allele)

    # Assign events to nodes
    # Each node must have at least one event (SNV or CNV)
    nodes = list(range(1,n_nodes)) + list(np.random.choice(range(1,n_nodes),n_SNVs+n_CNVs - n_nodes+1,replace=True))
    np.random.shuffle(nodes)
    
    # SNVs
    for i in range(n_SNVs):
        muts[nodes[i]].append(i)
    for n in range(n_nodes):
        muts[n] = sorted(muts[n])
    # CNVs
    loci_with_CNV = sorted(np.random.choice(n_SNVs,n_CNVs,replace=False))
    for i in range(n_CNVs):
        n = nodes[n_SNVs+i]
        loss_gain = np.random.choice([-1,1])
        allele = 0
        # copy number losses can only affect the reference allele,
        # and copy number gains can only affect the mutated alle when 
        # they are located in the same node as the corresponding mutation
        if loss_gain==1 and loci_with_CNV[i] in muts[n]:
            allele = np.random.choice([0,1])

        CNVs[n].append( (loci_with_CNV[i] , loss_gain , allele) )

    # Compute the genotypes of each node
    n_ref_allele = 2*np.ones((n_SNVs,n_nodes),dtype=int)
    n_alt_allele =  np.zeros((n_SNVs,n_nodes),dtype=int)
    for n in DFT_order:
        if n!=0:
            n_ref_allele[:,n] = n_ref_allele[:,parents[n]]
            n_alt_allele[:,n] = n_alt_allele[:,parents[n]]
        for mut in muts[n]:
            if n_ref_allele[mut,n]>0:
                n_ref_allele[mut,n]-=1
                n_alt_allele[mut,n]+=1
        for CNV in CNVs[n]:
            locus,gain_loss,allele = CNV
            if allele==0:
                n_ref_allele[locus,n]+=gain_loss
            else:
                if n_alt_allele[locus,n]==0:
                    raise Exception("Invalid CNV: affects the alt allele, but its copy number before CNV is 0")
                n_alt_allele[locus,n]+=gain_loss
            if (n_ref_allele[locus,n]<0 or n_alt_allele[locus,n]<0):
                raise Exception("Invalid CNV: copy number of one allele < 0.")


    # Sample node probabilities
    node_probabilities = np.random.dirichlet([nodeprobs_concentration]*n_nodes) 
    # Sample dropout rates
    dropout_rates = np.abs(np.random.normal(dropout_rate_avg,dropout_rate_sigma,n_SNVs))
    # Sample region probabilities
    region_probabilities = np.abs(np.random.normal(1,regionprob_sigma,n_SNVs))
    region_probabilities = [max(0.3,x) for x in region_probabilities]
    # Assign cells to nodes and generate read counts
    ref_reads =  np.zeros((n_SNVs,n_cells),dtype=int)
    alt_reads =  np.zeros((n_SNVs,n_cells),dtype=int)
    for j in range(n_cells):
        node = int(np.random.choice(n_nodes, p=node_probabilities))
        for i in range(n_SNVs):
            expected_depth = sequencing_depth * region_probabilities[i]*(n_ref_allele[i,node]+n_alt_allele[i,node])/2.0
            rate = np.random.gamma(shape=theta,scale=expected_depth/theta)
            depth = np.random.poisson(rate)
            c_r=0 # number of copies of the ref allele that did not get dropped out
            c_a=0 # number of copies of the alt allele that did not get dropped out
            for x in range(n_ref_allele[i,node]):
                if np.random.rand()>=dropout_rates[i]: c_r+=1
            for x in range(n_alt_allele[i,node]):
                if (np.random.rand())>= dropout_rates[i]: c_a+=1
            if (c_r==0 and c_a==0):
                if n_ref_allele[i,node]>0 and n_alt_allele[i,node]:
                    if np.random.choice(2)==0: c_r+=1
                    else: c_a+=1
                elif n_ref_allele[i,node]>0: c_r+=1
                else: c_a+=1
            f = c_a / (c_r+c_a) * (1-sequencing_error_rate) + c_r/(c_r+c_a) * sequencing_error_rate
            omega=1000 # high omega --> beta binomial ~= binomial
            a = f*omega
            b = (1-f)*omega
            alt_reads[i,j] = betabinom.rvs(depth,a,b)
            ref_reads[i,j] = depth-alt_reads[i,j]



    # Create files
    directory = os.path.dirname(basename)    
    if directory!="" and (not os.path.exists(directory)):
        os.makedirs(directory)

    d={}
    d["CHR"] = list(range(n_SNVs))
    d["REGION"] = ["Region"+str(x) for x in range(n_SNVs)]
    d["NAME"] = ["Mut"+str(x) for x in range(n_SNVs)]
    for j in range(n_cells):
        d[j]=[]
    for i in range(n_SNVs):
        for j in range(n_cells):
            d[j].append(str(ref_reads[i,j])+":"+str(alt_reads[i,j]))
    df_variants = pd.DataFrame(d)
    df_variants.to_csv(basename+"_variants.csv",index=False)

    depth=[]
    for i in range(n_SNVs):
        l=[]
        for j in range(n_cells):
            l.append(ref_reads[i,j]+alt_reads[i,j])
        depth.append(l)
    
    regions = np.array(depth)
    df_regions = pd.DataFrame(depth,index = [str(i)+"_Region"+str(i) for i in range(n_SNVs)]).astype(int)
    df_regions.to_csv(basename+"_regions.csv",index=True,header=False)
    #np.savetxt(basename+"_regions.csv",regions.astype(int),delimiter=",",fmt='%i')

    #Input for BiTSC2
    #np.savetxt(basename+"_DP.csv",regions.astype(int),delimiter=",",fmt='%i')
    #np.savetxt(basename+"_AD.csv",alt_reads.astype(int),delimiter=",",fmt='%i')

    # Save tree structure, with dropout rates and node probabilities
    parents_array = np.array(parents)
    np.savetxt(basename+"_tree.csv",parents_array,delimiter=",",fmt='%i')
    colors = ["lightcoral","skyblue3","sandybrown","paleturquoise3","thistle","darkolivegreen3","lightpink","mediumpurple","darkseagreen3","navajowhite","gold"]
    with open(basename+"_tree.gv", "w") as file:
        file.write("digraph G{\n")
        file.write("node [color=dimgray fontsize=24 fontcolor=black fontname=Helvetica penwidth=5];\n")
        for n in range(1,n_nodes):
            file.write(str(parents[n])+" -> "+str(n)+" [color=dimgray penwidth=4 weight=2];\n")
        for n in range(n_nodes):
            label = str(n)+"[label=<"
            for mut in muts[n]:
                label+=str(mut)+":Mut"+str(mut)+"<br/>"
            for CNV in CNVs[n]:
                label+="CNV"
                if CNV[1]>0: label+="+1"
                else: label+="-1"
                label+= " "+str(CNV[0])+":"+str(CNV[2])+"<br/>"

            if len(muts[n])==0 and len(CNVs[n])==0: label+=" "
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
    parser.add_argument('--nSNVs', type = int,default = 6, help='Number of SNVs')
    parser.add_argument('--nCNVs', type = int,default = 2, help='Number of CNVs')
    parser.add_argument('--dropoutsigma', type = float,default = 0.03, help='Standard deviation when sampling the dropout rates')
    parser.add_argument('--regionprobsigma', type = float,default = 0.00, help='Standard deviation when sampling the region probabilities')
    parser.add_argument('--theta', type = float,default = 10, help='Inverse overdispersion parameter for the negative binomial (when sampling the sequencing depth)')
    parser.add_argument('--depth', type = int,default = 10, help='Sequencing depth')
    parser.add_argument('--nodeprobconcentration', type=float, default=1,help='Concentration parameter of the dirichlet distribution when sampling the node probabilities.')
    args = parser.parse_args()
    generate_data(args.o,n_nodes=args.nnodes,n_cells=args.ncells,n_SNVs=args.nSNVs,n_CNVs=args.nCNVs,nodeprobs_concentration=args.nodeprobconcentration,\
        dropout_rate_avg=0.04,dropout_rate_sigma=args.dropoutsigma,seed=args.seed,sequencing_depth=args.depth,theta=args.theta,regionprob_sigma=args.regionprobsigma)