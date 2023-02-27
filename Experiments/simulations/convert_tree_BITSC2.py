import argparse
import os
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, help='Basename of the BITSC2 output')
parser.add_argument('-o', type = str, help='Output basename')
parser.add_argument('--metadata', type = str,help='Basename of the COMPASS/BITSC2 input (for metadata)')
parser.add_argument('--writeindices', action='store_true')
parser.set_defaults(writeindices=False)
args = parser.parse_args()

# Read parents
with open(args.i+"_treeBITSC2.csv") as infile:
    parents=[]
    for line in infile:
        parents.append(int(line)-1)
n_nodes = len(parents)


# Read metadata
df = pd.read_csv(args.metadata+"_variants.csv",sep=",")
n_variants = df.shape[0]
#SNV_to_region = [x.lstrip("Region") for x in df["REGION"]]
SNV_to_region = [x for x in df["REGION"]]
if "NAME" in df.columns:
    names = list(df["NAME"])
else:
    names = [str(x) for x in range(df.shape[0])]

df = pd.read_csv(args.metadata+"_regions.csv",sep=",",header=None)
region_to_index = {}
for i in range(df.shape[0]):
    region_to_index[df.iloc[i,0].split("_")[1]] = i

print(names)
print(SNV_to_region)
print(region_to_index)

# Add the regions without variants, which are encoded by variants after n_variants
regionsInd_with_variants = [region_to_index[region] for region in SNV_to_region]
for i in range(df.shape[0]):
    if not i in regionsInd_with_variants:
        SNV_to_region.append(df.iloc[i,0].split("_")[1])
        region_to_index[df.iloc[i,0].split("_")[1]] = i


# Reads SNVs and CNAs
SNVs = [[] for i in range(n_nodes)]
CNAs=[[] for i in range(n_nodes)]
SNV = 0
regions_with_CNA = set()
with open(args.i+"_originsBITSC2.csv","r") as file:
    for line in file:
        linesplit = line.split(" ")
        node = int(linesplit[0]) -1
        if SNV<n_variants:
            SNVs[node].append(SNV)
        nodeCNA= int(linesplit[2]) -1
        sign = np.sign(float(linesplit[3]))
        allele = int(linesplit[4])
        if sign!=0:
            print(SNV)
            print(len(SNV_to_region))
            region = SNV_to_region[SNV]
            if not region in regions_with_CNA: # each region can only contain one CNA (in case several loci are in the same region)
                CNAs[nodeCNA].append((region,sign,[]))
                regions_with_CNA.add(region)
            # Add the allele (there might be several variants in the region)
            for i in range(len(CNAs[nodeCNA])):
                if CNAs[nodeCNA][i][0]==region:
                    CNAs[nodeCNA][i][2].append(allele)
        SNV+=1
    n_muts=SNV
print(parents)
print(SNVs)
print(CNAs)

# Read cell attachments. The output cell attachments will depend on whether some nodes were removed....
nodes_counts = [0 for x in range(n_nodes)]
n_cells=0
with open(args.i+"_attachmentsBITSC2.csv","r") as infile:
    for line in infile:
        node = int(line.rstrip("\n"))-1
        nodes_counts[node]+=1
        n_cells+=1

# Remove empty nodes (and nodes with no cells attached )
node_mapping = [0 for x in nodes_counts] # map node in the original tree to nodes in the new tree (potentially with some nodes removed)
n=1 # node index in the new tree (possibly smaller)
m=1 # node index in the original tree
while n<n_nodes:
    if (len(SNVs[n])==0 and len(CNAs[n])==0): #or nodes_counts[n]==0:
        node_mapping[m] = parents[n]
        for m2 in range(len(node_mapping)):
            if node_mapping[m2]==n: node_mapping[m2]=parents[n]
        del SNVs[n]
        del CNAs[n]
        del nodes_counts[n]
        for n2 in range(1,n_nodes):
            if parents[n2]==n: parents[n2]=parents[n]
        del parents[n]
        n_nodes-=1
        # All nodes after n have to decrease their index by 1
        for n2 in range(1,n_nodes):
            if parents[n2]>n: parents[n2]-=1
        for m2 in range(len(node_mapping)):
            if node_mapping[m2]>n: node_mapping[m2]-=1
    else:
        node_mapping[m] = n
        n+=1
    m+=1
    
print("After removing empty nodes:")
print(parents)
print(SNVs)
print(CNAs)

# Write cell attachments.
nodes_counts_new = [0 for x in range(n_nodes)]
n_cells=0
with open(args.i+"_attachmentsBITSC2.csv","r") as infile:
    with open(args.o+"_cellAssignments.tsv","w") as outfile:
        outfile.write("cell\tnode\n")
        for line in infile:
            node = node_mapping[int(line.rstrip("\n"))-1]
            nodes_counts_new[node]+=1
            outfile.write(str(n_cells)+"\t"+str(node)+"\n")
            n_cells+=1


# Write
with open(args.o+"_tree.gv","w") as outfile:
    tmp = outfile.write("digraph G {\n")
    tmp = outfile.write("node [color=dimgray fontsize=24 fontcolor=black fontname=Helvetica penwidth=5];\n")

    # Write the parents
    for i in range(1,len(parents)):
        tmp = outfile.write(str(parents[i]) + " -> "+ str(i) + " [color=dimgray penwidth=4 weight=2];\n")

    # Write the content of each node
    tmp = outfile.write("0[label=< >];\n")
    for i in range(1,len(parents)):
        tmp = outfile.write(str(i)+"[label=<")
        for SNV in SNVs[i]:
            if args.writeindices:
                tmp = outfile.write(str(SNV)+": ")
            if SNV<len(names):
                tmp = outfile.write(names[SNV]+"<br/>")
            else:
                tmp = outfile.write(str(SNV)+"<br/>") # SNV corresponding to an empty region
        for CNA in CNAs[i]:
            if CNA[1]>0:
                tmp = outfile.write("Gain ")
            else:
                 tmp = outfile.write("Loss ")
            if args.writeindices:
                tmp = outfile.write(str(region_to_index[CNA[0]])+":")
            alleles = ["REF" if x==0 else "ALT" for x in CNA[2]]
            tmp = outfile.write(CNA[0]+":"+",".join(alleles)+"<br/>")
                
        outfile.write(">];\n")

    # Write the number of cells attached to each node
    colors = ["lightcoral","skyblue3","sandybrown","paleturquoise3","thistle","darkolivegreen3","lightpink","mediumpurple","darkseagreen3","navajowhite","gold"]
    for i in range(n_nodes):
        tmp = outfile.write(str(i) + " -> " + str(i+n_nodes)+ " [dir=none style=dashed weight=1 penwidth=5 color="+colors[i%len(colors)]+"];\n")
    
    for i in range(n_nodes):
        tmp = outfile.write(str(i+n_nodes)+"[label=\""+str(nodes_counts_new[i])+" cells\\n"+str(int(nodes_counts_new[i]*100//n_cells))+"%\" style = filled width=" +str(np.sqrt(100.0*nodes_counts[i]/n_cells) /3.0)+" color="+colors[i%len(colors)]+"];\n")

    outfile.write("}")