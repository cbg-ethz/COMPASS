# This script converts the true tree (in gv format), the tree inferred by COMPASS (gv) and the tree inferred by BITSC² to the format required
# by the tool that computes the MP3 distance.

import argparse
import os
import numpy as np
import pandas as pd
import copy

import mp3treesim as mp3

def map_SNV_regions(filename):
    SNV_to_region = []
    with open(filename,"r") as f:
        tmp = f.readline()
        for line in f:
            linesplit = line.split(",")
            region = int(linesplit[1].lstrip("Region"))
            SNV_to_region.append(region)
    return SNV_to_region


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

def read_tree_COMPASS(filename):
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
        tree["CNLOHs"] = [[] for x in tree["parents"]]
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
                        region = int(event[5:event.find(":")])+1
                        tree["CNLOHs"][node].append(region)
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
    SNV_to_region = map_SNV_regions(basename+"_variants.csv")
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
    regions_with_CNV = set()
    with open(basename+"_originsBITSC2.csv","r") as file:
        for line in file:
            linesplit = line.split(" ")
            node = int(linesplit[0]) -1
            tree["muts"][node].append(mut)
            nodeCNV= int(linesplit[2]) -1
            sign = np.sign(float(linesplit[3]))
            if sign!=0:
                region = SNV_to_region[mut-1]+1
                if not region in regions_with_CNV: # each region can only contain one CNV
                    tree["CNVs"][nodeCNV].append((region,sign))
                    regions_with_CNV.add(region)
            mut+=1
        tree["n_muts"]=mut
    remove_empty_nodes(tree)
    return tree

def read_tree_infSCITE(filename):
    """ Read a tree stored as a graphviz file, output by SCITE."""
    with open(filename,"r") as file:
        tree={}
        tree["parents"]=[]
        # skip first lines
        tmp = file.readline()
        tmp = file.readline()
        # Skip node labels
        finished_reading_labels=False
        while not finished_reading_labels:
            line = file.readline()
            if not line.find("[")>0:
                finished_reading_labels = True
        # Read parents
        finished_reading_parents=False
        while not finished_reading_parents:
            line_split = line.rstrip(";\n").split("->")
            parent = int(line_split[0])
            child = int(line_split[1])
            if max(parent,child)+1 > len(tree["parents"]):
                tree["parents"]+=[-1]* (max(parent,child)+ 1 - len(tree["parents"]))
            tree["parents"][child] = parent
            line = file.readline()
            if not line.find("->")>0:
                finished_reading_parents=True

        #set root to 0
        tree["parents"] = [tree["parents"][-1]] + tree["parents"][:-1]
        for i in range(1,len(tree["parents"])):
            if tree["parents"][i]== len(tree["parents"])-1:
                tree["parents"][i]=0
            else:
                tree["parents"][i]+=1
        tree["muts"] = [[]]
        tree["CNVs"] = [[] for x in tree["parents"]]
        tree["CNLOHs"] = [[] for x in tree["parents"]]
        tree["n_muts"]=0
        for i in range(1,len(tree["parents"])):
            tree["muts"].append([i])
            tree["n_muts"]+=1
    remove_empty_nodes(tree)
    return tree




def get_node_label_MP3(tree,node,CNVs_map):
    label_list=[]
    for mut in tree["muts"][node]:
        label_list.append(str(mut))
    for CNV in tree["CNVs"][node]:
        if not CNV in CNVs_map:
            CNVs_map[CNV] = tree["n_muts"]+1+len(CNVs_map)
        label_list.append(str(CNVs_map[CNV]))
    return ",".join(label_list)



def write_tree_MP3(tree,CNVs_map,filename):
    with open(filename,"w") as file:
        file.write("digraph Tree { \n")
        for n in range(len(tree["parents"])):
            file.write(str(n) + " [label=\""+get_node_label_MP3(tree,n,CNVs_map)+"\"];\n")
        for n in range(len(tree["parents"])):
            file.write(str(tree["parents"][n]) + " -> " + str(n) +";\n")
        file.write("}")


def convert_trees(basename_in, basename_tmp,method):
    tree_true = read_tree_COMPASS(basename_in+"_treeTRUE.gv")
    if method =="COMPASS":
        tree_inferred = read_tree_COMPASS(basename_in+"_treeCOMPASS.gv")
    elif method=="SCITE":
        tree_inferred = read_tree_infSCITE(basename_in+"_map0.gv")
    else:
        tree_inferred = read_tree_BITSC2(basename_in)
    CNVs_map = {}
    write_tree_MP3(tree_true,CNVs_map,basename_tmp+"TRUE-MP3.gv")
    write_tree_MP3(tree_inferred,CNVs_map,basename_tmp+method+"-MP3.gv")

    
def evaluate_distances(basename_tmp,output,method):
    treeTRUE = mp3.read_dotfile(basename_tmp+"TRUE-MP3.gv")
    treeINFERRED = mp3.read_dotfile(basename_tmp+method+"-MP3.gv")
    with open(output,"w") as outfile:
        outfile.write(str(mp3.similarity(treeTRUE,treeINFERRED)))


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type = str, help='Input basename')
    parser.add_argument('-m', type = str, help='Name of the method (COMPASS, BITSC2 or SCITE')
    parser.add_argument('--tmp', type = str, help='Temp files basename')
    parser.add_argument('-o', type = str, help='Output file')
    args = parser.parse_args()
    convert_trees(args.i,args.tmp,args.m)
    evaluate_distances(args.tmp,args.o,args.m)
    