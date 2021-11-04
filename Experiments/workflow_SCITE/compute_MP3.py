# This script converts the true tree (in gv format), the tree inferred by COMPASS (gv) and the tree inferred by BITSC² to the format required
# by the tool that computes the Bourque distance.

import argparse
import os
import numpy as np
import pandas as pd
import copy

import mp3treesim as mp3

def remove_empty_nodes(tree):
    if len(tree["muts"][0])==0:
        tree["muts"][0].append(0)
    node = 1
    size=len(tree["parents"])
    while node < size:
        if len(tree["muts"][node])==0:
            del tree["muts"][node]
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
                        pass

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

    remove_empty_nodes(tree)
    return tree

def read_tree_gv_SCITE(filename):
    """ Read a tree stored as a graphviz file, output by SCITE."""
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
                line_split = line.rstrip(";\n").split("->")
                parent = int(line_split[0])-1
                child = int(line_split[1])-1
                if max(parent,child)+1 > len(tree["parents"]):
                    tree["parents"]+=[-1]* (max(parent,child)+ 1 - len(tree["parents"]))
                tree["parents"][child] = parent
            else:
                finished_reading_parents=True
        #set root to 0
        tree["parents"] = [tree["parents"][-1]] + tree["parents"][:-1]
        for i in range(1,len(tree["parents"])):
            if tree["parents"][i]== len(tree["parents"])-1:
                tree["parents"][i]=0
            else:
                tree["parents"][i]+=1
        tree["muts"] = [[]]
        for i in range(1,len(tree["parents"])):
            tree["muts"].append([i])
    remove_empty_nodes(tree)
    return tree



def get_node_label(tree,node):
    label_list=[]
    for mut in tree["muts"][node]:
        label_list.append(str(mut))
    return "_".join(label_list)

def get_node_label_MP3(tree,node):
    label_list=[]
    for mut in tree["muts"][node]:
        label_list.append(str(mut))
    return ",".join(label_list)

def write_tree(tree,file):
    for n in range(1,len(tree["parents"])):
        if tree["parents"][n]==0:
            label = get_node_label(tree,n)
            label_parent = get_node_label(tree,tree["parents"][n])
            file.write(label_parent+" " + label + "\n")
    for n in range(1,len(tree["parents"])):
        if tree["parents"][n]!=0:
            label = get_node_label(tree,n)
            label_parent = get_node_label(tree,tree["parents"][n])
            file.write(label_parent+" " + label + "\n")

def write_tree_MP3(tree,filename):
    with open(filename,"w") as file:
        file.write("digraph Tree { \n")
        for n in range(len(tree["parents"])):
            file.write(str(n) + " [label=\""+get_node_label_MP3(tree,n)+"\"];\n")
        for n in range(len(tree["parents"])):
            file.write(str(tree["parents"][n]) + " -> " + str(n) +";")
        file.write("}")


def convert_trees(basename_in, basename_tmp):
    tree_true = read_tree_gv(basename_in+"_treeTRUE.gv")
    tree_COMPASS = read_tree_gv(basename_in+"_treeCOMPASS.gv")
    tree_SCITE = read_tree_gv_SCITE(basename_in+"_map0.gv")
    write_tree_MP3(tree_true,basename_tmp+"TRUE-MP3.gv")
    write_tree_MP3(tree_COMPASS,basename_tmp+"COMPASS-MP3.gv")
    write_tree_MP3(tree_SCITE,basename_tmp+"SCITE-MP3.gv")

    
def evaluate_distances(basename_tmp,output):
    treeTRUE = mp3.read_dotfile(basename_tmp+"TRUE-MP3.gv")
    treeCOMPASS = mp3.read_dotfile(basename_tmp+"COMPASS-MP3.gv")
    treeSCITE = mp3.read_dotfile(basename_tmp+"SCITE-MP3.gv")
    with open(output,"w") as outfile:
        outfile.write(str(mp3.similarity(treeTRUE,treeCOMPASS))+"\n")
        outfile.write(str(mp3.similarity(treeTRUE,treeSCITE)))

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type = str, help='Input basename')
    parser.add_argument('--tmp', type = str, help='Temp files basename')
    parser.add_argument('-o', type = str, help='Output file')
    args = parser.parse_args()
    convert_trees(args.i,args.tmp)
    evaluate_distances(args.tmp,args.o)
    