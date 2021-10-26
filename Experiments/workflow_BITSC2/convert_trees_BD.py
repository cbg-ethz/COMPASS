# This script converts the true tree (in gv format), the tree inferred by COMPASS (gv) and the tree inferred by BITSC² to the format required
# by the tool that computes the Bourque distance.

import argparse
import os
import numpy as np
import pandas as pd
import copy

def remove_empty_nodes(tree):
    if len(tree["muts"][0])==0 and len(tree["CNVs"][0])==0:
        tree["muts"][0].append(0)
    node = 1
    size=len(tree["parents"])
    while node < size:
        if len(tree["muts"][node])==0 and len(tree["CNVs"][node])==0:
            del tree["muts"][node]
            del tree["CNVs"][node]
            parent = tree["parents"][node]
            del tree["parents"][node]
            if parent>node: parent-=1
            for i in range(len(tree["parents"])):
                if tree["parents"][i]==node: tree["parents"][i] = parent
                elif tree["parents"][i]>node: tree["parents"][i]-=1
            node=0
            size-=1
        else:
            node+=1

def read_tree_gv(filename):
    """ Read a tree stored as a graphviz file."""
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
        tree["muts"] = [[] for x in tree["parents"]]
        tree["CNVs"] = [[] for x in tree["parents"]]
        tree["n_muts"]=0
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
                    if event[:3]=="CNV":
                        sign = 1 if event[3]=="+" else -1
                        region = int(event[5:event.find(":")])+1
                        tree["CNVs"][node].append((region,sign))

                    elif event[:3]=="CNL":
                        pass
                    else:
                        end = event.find(":")
                        if end>0:
                            tree["muts"][node].append(int(event[:end])+1)
                            tree["n_muts"]+=1
                        else:
                            pass
                line = file.readline()

    # make sure that each event is present in only one node
    events_set = set()
    for n in range(len(tree["parents"])):
        l = copy.deepcopy(tree["CNVs"][n])
        for CNV in l:
            if CNV in events_set:
                tree["CNVs"][n].remove(CNV)
            else:
                events_set.add(CNV)
    remove_empty_nodes(tree)
    return tree

def read_tree_BITSC2(basename):
    tree={}
    with open(basename+"_treeBITSC2.csv","r") as file:
        parents=[]
        for line in file:
            parents.append(int(line)-1)
        tree["parents"] = parents
    n_nodes = len(tree["parents"])

    tree["muts"] = [[] for i in range(n_nodes)]
    tree["CNVs"]=[[] for i in range(n_nodes)]
    mut = 1
    with open(basename+"_originsBITSC2.csv","r") as file:
        for line in file:
            linesplit = line.split(" ")
            node = int(linesplit[0]) -1
            tree["muts"][node].append(mut)
            nodeCNV= int(linesplit[2]) -1
            sign = np.sign(float(linesplit[3]))
            if sign!=0:
                tree["CNVs"][nodeCNV].append((mut,sign))
            mut+=1
        tree["n_muts"]=mut
    remove_empty_nodes(tree)
    return tree



def get_node_label(tree,node,CNVs_map):
    label_list=[]
    for mut in tree["muts"][node]:
        label_list.append(str(mut))
    for CNV in tree["CNVs"][node]:
        if not CNV in CNVs_map:
            CNVs_map[CNV] = tree["n_muts"]+1+len(CNVs_map)
        label_list.append(str(CNVs_map[CNV]))
    return "_".join(label_list)

def write_tree(tree,CNVs_map,file):
    for n in range(1,len(tree["parents"])):
        if tree["parents"][n]==0:
            label = get_node_label(tree,n,CNVs_map)
            label_parent = get_node_label(tree,tree["parents"][n],CNVs_map)
            file.write(label_parent+" " + label + "\n")
    for n in range(1,len(tree["parents"])):
        if tree["parents"][n]!=0:
            label = get_node_label(tree,n,CNVs_map)
            label_parent = get_node_label(tree,tree["parents"][n],CNVs_map)
            file.write(label_parent+" " + label + "\n")


def convert_trees(basename_in, output):
    tree_true = read_tree_gv(basename_in+"_treeTRUE.gv")
    tree_COMPASS = read_tree_gv(basename_in+"_treeCOMPASS.gv")
    tree_BITSC2 = read_tree_BITSC2(basename_in)
    CNVs_map = {}
    with open(output,"w") as file:
        file.write("#tree TRUE: " + str(len(tree_true["muts"])-1) + "\n")
        write_tree(tree_true,CNVs_map,file)
        file.write("#tree COMPASS: " + str(len(tree_COMPASS["muts"])-1) + "\n")
        write_tree(tree_COMPASS,CNVs_map,file)
        file.write("#tree BITSC2: " + str(len(tree_BITSC2["muts"])-1) + "\n")
        write_tree(tree_BITSC2,CNVs_map,file)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type = str, help='Input basename')
    parser.add_argument('-o', type = str, help='Output file')
    args = parser.parse_args()
    convert_trees(args.i,args.o)
    