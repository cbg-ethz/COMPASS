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
    if len(tree["SNVs"][0])==0 and len(tree["CNAs"][0])==0:
        tree["SNVs"][0].append(0)
    node = 1
    size=len(tree["parents"])
    while node < size:
        if len(tree["SNVs"][node])==0 and len(tree["CNAs"][node])==0:
            del tree["SNVs"][node]
            del tree["CNAs"][node]
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

def read_tree(filename,use_SNV,use_CNA):
    """ Read a tree stored as a graphviz file. Optionally ignore SNVs or CNAs."""
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
        if not use_SNV:
            tree["SNVs"] = [[] for x in tree["parents"]]
            tree["nSNVs"]=0
        if not use_CNA:
            tree["CNAs"] = [[] for x in tree["parents"]]


    # make sure that each event is present in only one node
    events_set = set()
    for n in range(len(tree["parents"])):
        l = copy.deepcopy(tree["CNAs"][n])
        for CNA in l:
            if CNA in events_set:
                tree["CNAs"][n].remove(CNA)
            else:
                events_set.add(CNA)
    remove_empty_nodes(tree)
    return tree



def get_node_label_MP3(tree,node,CNAs_map):
    label_list=[]
    for SNV in tree["SNVs"][node]:
        label_list.append(str(SNV))
    for CNA in tree["CNAs"][node]:
        if not CNA in CNAs_map:
            CNAs_map[CNA] = tree["nSNVs"]+1+len(CNAs_map)
        label_list.append(str(CNAs_map[CNA]))
    return ",".join(label_list)



def write_tree_MP3(tree,CNAs_map,filename):
    with open(filename,"w") as file:
        file.write("digraph Tree { \n")
        for n in range(len(tree["parents"])):
            file.write(str(n) + " [label=\""+get_node_label_MP3(tree,n,CNAs_map)+"\"];\n")
        for n in range(1,len(tree["parents"])):
            file.write(str(tree["parents"][n]) + " -> " + str(n) +";\n")
        file.write("}")


def convert_trees(true_tree_file,inferred_tree_file, basename_tmp,useSNV,useCNA):
    tree_true = read_tree(true_tree_file,use_SNV=useSNV,use_CNA=useCNA)
    tree_inferred = read_tree(inferred_tree_file,use_SNV=useSNV,use_CNA=useCNA)
    CNAs_map = {}
    write_tree_MP3(tree_true,CNAs_map,basename_tmp+"TRUE-MP3.gv")
    write_tree_MP3(tree_inferred,CNAs_map,basename_tmp+"INFERRED-MP3.gv")

    
def evaluate_distances(basename_tmp,output,union):
    treeTRUE = mp3.read_dotfile(basename_tmp+"TRUE-MP3.gv")
    treeINFERRED = mp3.read_dotfile(basename_tmp+"INFERRED-MP3.gv")
    print("== MP3 ==")
    print(treeTRUE.label_set)
    print(treeTRUE.node_to_labels)
    print(treeINFERRED.label_set)
    print(treeINFERRED.node_to_labels)

    
    if union: # for CNAs, use MP3 in union mode, because with the sigmoid mode there are divisions by 0 in some cases.
        if len(treeINFERRED.label_set)==1 and len(treeTRUE.label_set)==1:
            print("Both trees are empty --> similarity = 1")
            with open(output,"w") as outfile:
                outfile.write("1.0")
        elif len(treeINFERRED.label_set)==1 or len(treeTRUE.label_set)==1:
            print("One tree is empty but the other is not --> similarity = 0")
            with open(output,"w") as outfile:
                outfile.write("0.0")
        else:
            with open(output,"w") as outfile:
                outfile.write(str(mp3.similarity(treeTRUE,treeINFERRED,mode="union")))
    else:
        if len(treeINFERRED.label_set)==1: # empty tree
           with open(output,"w") as outfile:
                outfile.write("0.0") 
        else:
            with open(output,"w") as outfile:
                outfile.write(str(mp3.similarity(treeTRUE,treeINFERRED,mode="sigmoid")))


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', type = str, help='True tree')
    parser.add_argument('-i2', type = str, help='Inferred tree')
    parser.add_argument('--tmp', type = str, help='Temp files basename')
    parser.add_argument('--useSNV', type = int,default=1, help='Consider SNVs in the MP3 similarity')
    parser.add_argument('--useCNA', type = int,default=1, help='Consider CNAs in the MP3 similarity')
    parser.add_argument('-o', type = str, help='Output file')
    args = parser.parse_args()
    useSNV = args.useSNV==1
    useCNA = args.useCNA==1
    print((args.o,useSNV,useCNA))
    convert_trees(args.i1,args.i2,args.tmp,useSNV,useCNA)
    evaluate_distances(args.tmp,args.o, not useSNV)
    