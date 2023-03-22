

import argparse
import os
import numpy as np
import pandas as pd
import copy


def read_tree(filename):
    """ Read a tree stored as a graphviz file."""
    with open(filename,"r") as infile: # Check that the tree is not empty.
        if len(infile.readlines())<=1:
            return {"parents":[-1],"SNVs":[[]],"CNAs":[[]],'n_SNVs':0}
    with open(filename,"r") as file:
        tree={}
        tree["parents"]=[]
        # skip first lines
        tmp = file.readline()
        tmp = file.readline()
        # Read parents
        finished_reading_parents=False
        while not finished_reading_parents:
            line = file.readline()
            if line.find("->")>0:
                line_split = (line.split("[")[0]).split("->")
                parent = int(line_split[0])
                child = int(line_split[1])
                if max(parent,child)+1 > len(tree["parents"]):
                    tree["parents"]+=[-1]* (max(parent,child)+ 1 - len(tree["parents"]))
                tree["parents"][child] = parent
            else:
                finished_reading_parents=True
        if tree["parents"]==[]: tree["parents"] = [-1]
        
        # Read node labels
        tree["SNVs"] = [[] for x in tree["parents"]]
        tree["CNAs"] = [[] for x in tree["parents"]]
        tree["nSNVs"]=0
        finished_reading_labels=False
        while not finished_reading_labels:
            if line.find("->")>=0:
                finished_reading_labels = True
            else:
                node = int(line[0:line.find("[")])
                linesplit=line[line.find("[")+8:line.find("]")-1].split("<br/>")
                if len(linesplit)>0: linesplit=linesplit[:-1]

                for event in linesplit:
                    event = event.lstrip("<B>").rstrip("</B>")
                    if event[:4]=="Gain":
                        region = int(event[5:event.find(":")])+1
                        tree["CNAs"][node].append((region,1))
                    elif event[:4]=="Loss":
                        region = int(event[5:event.find(":")])+1
                        tree["CNAs"][node].append((region,-1))
                    elif event[:3]=="CNL":
                        region = int(event[6:event.find(":")])+1
                        tree["CNAs"][node].append((region,0))
                    else:
                        end = event.find(":")
                        if end>0:
                            tree["SNVs"][node].append(int(event[:end])+1)
                            tree["nSNVs"]+=1
                        else:
                            pass
                line = file.readline()

    # make sure that each event is present in only one node
    events_set = set()
    for n in range(len(tree["parents"])):
        l = copy.deepcopy(tree["CNAs"][n])
        for CNA in l:
            if CNA in events_set:
                tree["CNAs"][n].remove(CNA)
            else:
                events_set.add(CNA)
    # Compute genotypes of each node
    children = [[] for x in range(len(tree["parents"]))]
    for i in range(1,len(tree["parents"])):
        children[tree["parents"][i]].append(i)
    DFT_order=[]
    stack = [0]
    while stack!=[]:
        top = stack.pop()
        DFT_order.append(top)
        for child in children[top]: stack.append(child)
    

    genotypes = [set() for x in range(len(tree["parents"]))]
    for n in DFT_order:
        if tree["parents"][n]>=0:
            gen = copy.copy(genotypes[tree["parents"][n]])
            for x in tree["SNVs"][n]: gen.add(x)
            for x in tree["CNAs"][n]: gen.add(x)
            genotypes[n] = gen
    tree["genotypes"] = genotypes

    return tree


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', type = str, help='True tree')
    parser.add_argument('-i2', type = str, help='Inferred tree')
    parser.add_argument('-a1', type = str, help='True cell assignments')
    parser.add_argument('-a2', type = str, help='Inferred cell assignments')
    parser.add_argument('-o', type = str, help='Output file')
    args = parser.parse_args()
    print("-----------------------------------")
    print("Computing cell assignment accuracy")
    tree_true=read_tree(args.i1)
    tree_inferred=read_tree(args.i2)
    if len(tree_inferred["parents"])==1: # only one node: consider that all cells are mapped to the root
        trueAssignments = pd.read_csv(args.a1,sep="\t",index_col=0)
        accuracy = np.sum(trueAssignments["node"]==0) / trueAssignments.shape[0]
        print("Empty tree; using as accuracy the number of cells assigned to the root in the true tree:" + str(accuracy))
        with open(args.o,"w") as outfile:
            outfile.write(str(accuracy))
    else:
        print(tree_true["genotypes"])
        print(tree_inferred["genotypes"])
        # Map each node of the inferred tree to the node with the closest genotype in the true tree
        node_mapping=[]
        for n in range(len(tree_inferred["parents"])):
            min_dist=10000000
            best_node=0
            for m in range(len(tree_true["parents"])):
                dist = len(tree_inferred["genotypes"][n].difference(tree_true["genotypes"][m])) + len(tree_true["genotypes"][m].difference(tree_inferred["genotypes"][n]))
                if dist < min_dist:
                    min_dist = dist
                    best_node = m
            node_mapping.append(best_node)
        print(node_mapping)


        trueAssignments = pd.read_csv(args.a1,sep="\t",index_col=0)
        inferredAssignments = pd.read_csv(args.a2,sep="\t",index_col=0)
        n_matches=0
        for i in inferredAssignments.index:
            if node_mapping[inferredAssignments.loc[i,"node"]]==trueAssignments.loc[i,"node"]:
                n_matches+=1
        accuracy = n_matches / inferredAssignments.shape[0]
        with open(args.o,"w") as outfile:
                outfile.write(str(accuracy))

    
    