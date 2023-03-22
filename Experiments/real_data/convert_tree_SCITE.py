import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, help='SCITE tree in .gv format')
parser.add_argument('-o', type = str, help='Output basename')
parser.add_argument('--names', type = str,default = None, help='Mutation names')
parser.add_argument('--writeindices', action='store_true')
parser.add_argument('--removeemptynodes', action='store_true', help='If set, remove nodes with few cells attached.')
parser.set_defaults(writeindices=False)
parser.set_defaults(removeemptynodes=False)
args = parser.parse_args()



with open(args.i,"r") as infile:
    with open(args.o+"_tree.gv","w") as outfile:
        tmp = outfile.write(infile.readline())
        tmp = infile.readline()
        tmp = outfile.write("node [color=dimgray fontsize=24 fontcolor=black fontname=Helvetica penwidth=5];\n")
        
        # Read the tree
        finished_reading_parents=False
        parents=[]
        # Read parents. SCITE starts nodes at 1; infSCITE at 0
        while not finished_reading_parents:
            line = infile.readline()
            if line.find("->")>0:
                line_split = line.rstrip(";\n").split("->")
                parent = int(line_split[0])
                child = int(line_split[1])
                if max(parent,child)+1 > len(parents):
                    parents+=[-1]* (max(parent,child)+ 1 - len(parents))
                parents[child] = parent
            else:
                finished_reading_parents=True

        # Read node attachments
        node_cellcounts=[0 for x in parents]
        cell2node={}
        finished_reading_attachments = False
        while not finished_reading_attachments:
            line = infile.readline()
            if line.find("->")>0:
                line_split = line.rstrip(";\n").split("->")
                node = int(line_split[0])+1
                if node==len(parents): node=0
                cell = line_split[1]
                if cell[0]==" ": cell=cell[1:]
                if cell[0]=="s": cell=cell[1:]
                if cell[0]=="_": cell = cell[1:]
                cell = int(cell)
                cell2node[cell] = node
                node_cellcounts[node]+=1
            else:
                finished_reading_attachments=True
        n_cells = 0
        for x in node_cellcounts:
            n_cells+=x
        
        print(node_cellcounts)

        #set root to 0
        parents = [parents[-1]] + parents[:-1]
        for i in range(1,len(parents)):
            if parents[i]== len(parents)-1:
                parents[i]=0
            else:
                parents[i]+=1
        SNVs = [[]]
        for i in range(1,len(parents)):
            SNVs.append([i-1])
        
        n_nodes_initial = len(parents)

        

        # Remove empty with few cells attached to them, which have exactly one child, and move the SNVs to the child.
        if args.removeemptynodes:
            children=[[] for x in range(len(parents))]
            for n in range(1,len(parents)):
                children[parents[n]].append(n)

            
            node_mapping = [0 for x in range(len(parents))] # map node in the original tree to nodes in the new tree (potentially with some nodes removed)
            n=1 # node index in the new tree (possibly smaller)
            m=1 # node index in the original tree
            while n < len(parents):
                if node_cellcounts[n] <= 0.02*n_cells and len(children[n])==1:
                    SNVs[children[n][0]] += SNVs[n]
                    node_mapping[m] = children[n][0]
                    for m2 in range(1,len(node_mapping)):
                        if node_mapping[m2]==n: node_mapping[m2] = children[n][0]
                    del SNVs[n]
                    del node_cellcounts[n]
                    for n2 in range(1,len(parents)):
                        if parents[n2]==n: parents[n2]=parents[n]
                    for n2 in range(1,len(parents)):
                        if parents[n2]>n: parents[n2]-=1
                    del parents[n]
                    children=[[] for x in range(len(parents))]
                    for n2 in range(1,len(parents)):
                        children[parents[n2]].append(n2)
                    for m2 in range(len(node_mapping)):
                        if node_mapping[m2]>n: node_mapping[m2]-=1
                else:
                    node_mapping[m]=n
                    n+=1
                m+=1
        else:
            node_mapping = [x for x in range(len(parents))]
        print(node_mapping)
        SNVs = [sorted(l) for l in SNVs]

        # write cell assignments
        cells=[]
        nodes=[]
        nodes_counts_new = [0 for x in parents]
        print(parents)
        print(nodes_counts_new)
        for x in sorted(cell2node.keys()):
            cells.append(x)
            #print((x,cell2node[x],node_mapping[cell2node[x]]))
            nodes.append(node_mapping[cell2node[x]])
            nodes_counts_new[node_mapping[cell2node[x]]]+=1
        df_assignments = pd.DataFrame({"cell":cells,"node":nodes})
        df_assignments.to_csv(args.o+"_cellAssignments.tsv",sep="\t",index=False)


        # Write the parents
        for i in range(1,len(parents)):
            tmp = outfile.write(str(parents[i]) + " -> "+ str(i) + " [color=dimgray penwidth=4 weight=2];\n")

        if args.names is not None:
            names=[]
            with open(args.names,"r") as infile:
                for line in infile:
                    names.append(line.rstrip("\n"))
        else:
            names = ["SNV"+str(x) for x in range(n_nodes_initial-1)]
        print(names)
        print(SNVs)

        # Write the content of each node
        tmp = outfile.write("0[label=< >];\n")
        for i in range(1,len(parents)):
            tmp = outfile.write(str(i)+"[label=<")
            for SNV in SNVs[i]:
                print((i,SNV))
                if args.writeindices:
                    tmp = outfile.write(str(SNV)+": "+names[SNV]+"<br/>")
                else:
                    tmp = outfile.write(names[SNV]+"<br/>")
            outfile.write(">];\n")

        # Write the cell attachments counts
        n_nodes = len(parents)
        n_cells = 0
        for x in nodes_counts_new:
            n_cells+=x
        colors = ["lightcoral","skyblue3","sandybrown","paleturquoise3","thistle","darkolivegreen3","lightpink","mediumpurple","darkseagreen3","navajowhite","gold"]
        for i in range(n_nodes):
            tmp = outfile.write(str(i) + " -> " + str(i+n_nodes)+ " [dir=none style=dashed weight=1 penwidth=5 color="+colors[i%len(colors)]+"];\n")
        
        for i in range(n_nodes):
            tmp = outfile.write(str(i+n_nodes)+"[label=\""+str(nodes_counts_new[i])+" cells\\n"+str(int(nodes_counts_new[i]*100//n_cells))+"%\" style = filled width=" +str(np.sqrt(100.0*nodes_counts_new[i]/n_cells) /3.0)+" color="+colors[i%len(colors)]+"];\n")
        outfile.write("}")

