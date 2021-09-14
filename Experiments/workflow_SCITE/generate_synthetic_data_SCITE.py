import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import nbinom, betabinom


def generate_data(basename,n_muts,n_cells,nodeprobs_concentration,dropout_rate_avg,dropout_rate_sigma,seed=0):
    np.random.seed(seed)
    n_nodes = n_muts+1
    sequencing_error_rate = 0.01
    omega_hom = 20
    omega_het = 6
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


    # Add mutations to the nodes
    muts = [[] for x in range(n_nodes)]
    for i in range(1,n_nodes): muts[i] = [i-1]


    
    # Compute the genotypes of each node
    n_ref_allele = 2*np.ones((n_muts,n_nodes),dtype=int)
    n_alt_allele =  np.zeros((n_muts,n_nodes),dtype=int)
    for n in DFT_order:
        if n!=0:
            n_ref_allele[:,n] = n_ref_allele[:,parents[n]]
            n_alt_allele[:,n] = n_alt_allele[:,parents[n]]
        for mut in muts[n]:
            if n_ref_allele[mut,n]>0:
                n_ref_allele[mut,n]-=1
                n_alt_allele[mut,n]+=1

    # Sample node probabilities (minimum node probability of 2%)
    node_probabilities = np.random.dirichlet([nodeprobs_concentration]*n_nodes) +0.02
    node_probabilities = node_probabilities / sum(node_probabilities)
    # Sample dropout rates
    dropout_rates = np.abs(np.random.normal(dropout_rate_avg,dropout_rate_sigma,n_muts))
    # Assign cells to nodes and generate read counts
    ref_reads =  np.zeros((n_muts,n_cells),dtype=int)
    alt_reads =  np.zeros((n_muts,n_cells),dtype=int)
    genotypes = np.zeros((n_muts,n_cells),dtype=int)
    for j in range(n_cells):
        node = int(np.random.choice(n_nodes, p=node_probabilities))
        for i in range(n_muts):
            depth = np.random.randint(10,40)
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
            if c_r==0 or c_a==0: omega = omega_hom
            else: omega = omega_het
            a = f*omega
            b = (1-f)*omega
            alt_reads[i,j] = betabinom.rvs(depth,a,b)
            ref_reads[i,j] = depth-alt_reads[i,j]


            if ref_reads[i,j] + alt_reads[i,j]<8:
                genotypes[i,j]=3
            elif min(ref_reads[i,j],alt_reads[i,j])>0.17 * (ref_reads[i,j] + alt_reads[i,j]): #15
                genotypes[i,j]=1
            elif alt_reads[i,j]>0.88 * (alt_reads[i,j] + ref_reads[i,j]): #90
                genotypes[i,j]=2
            elif ref_reads[i,j]>0.88 * (alt_reads[i,j] + ref_reads[i,j]):
                genotypes[i,j]=0
            else:
                genotypes[i,j]=3

    # Create files
    directory = os.path.dirname(basename)    
    if directory!="" and (not os.path.exists(directory)):
        os.makedirs(directory)

    np.savetxt(basename+"_genotypes.csv",genotypes.astype(int),delimiter=" ",fmt='%i')

    d={}
    d["CHR"] = list(range(n_muts))
    d["REGION"] = list(range(n_muts))
    for j in range(n_cells):
        d[j]=[]
    for i in range(n_muts):
        for j in range(n_cells):
            d[j].append(str(ref_reads[i,j])+":"+str(alt_reads[i,j]))
    df_variants = pd.DataFrame(d)
    df_variants.to_csv(basename+"_variants.csv",index=False)


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
                label+=str(mut)+":"+str(mut)+"<br/>"
            if len(muts[n])==0: label+=" "
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
    parser.add_argument('--ncells', type = int,default = 2000, help='Number of cells')
    parser.add_argument('--nmuts', type = int,default = 6, help='Number of mutations')
    parser.add_argument('--dropoutrate', type = float,default = 0.04, help='Dropout rate')
    parser.add_argument('--dropoutsigma', type = float,default = 0.0, help='Standard deviation when sampling the dropout rates')
    parser.add_argument('--nodeprobconcentration', type=float, default=1,help='Concentration parameter of the dirichlet distribution when sampling the node probabilities.')
    args = parser.parse_args()
    generate_data(args.o,n_cells=args.ncells,n_muts=args.nmuts,nodeprobs_concentration=args.nodeprobconcentration,\
        dropout_rate_avg=args.dropoutrate,dropout_rate_sigma=args.dropoutsigma,seed=args.seed)