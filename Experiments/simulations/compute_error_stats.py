import argparse
import os
import numpy as np
import pandas as pd
import copy

SNV_to_region = {}

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

def read_tree(filename,use_SNV=True,use_CNA=True):
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



def count_FP_TN(treeTRUE,treeINFERRED):
    true_CNAs=[]
    for CNAs in treeTRUE["CNAs"]:
        for CNA in CNAs:
            true_CNAs.append(CNA)
    inferred_CNAs=[]
    for CNAs in treeINFERRED["CNAs"]:
        for CNA in CNAs:
            inferred_CNAs.append(CNA)
    false_positives = 0
    true_negatives = 0
    for CNA in inferred_CNAs:
        if not CNA in true_CNAs:
            false_positives+=1
    n_regions = 0
    for SNVs in treeTRUE["SNVs"]:
        for SNV in SNVs:
            n_regions = max(n_regions,SNV+1)
    for i in range(n_regions):
        is_negative=True
        for CNA in true_CNAs:
            if CNA[0]==i:
                is_negative=False
        if is_negative:
            not_inferred=True
            for CNA in inferred_CNAs:
                if CNA[0]==i:
                    not_inferred=False
            if not_inferred:
                true_negatives+=1
    return (false_positives,true_negatives)

def count_TP_FN(treeTRUE,treeINFERRED):
    true_CNAs=[]
    for CNAs in treeTRUE["CNAs"]:
        for CNA in CNAs:
            true_CNAs.append(CNA)
    inferred_CNAs=[]
    for CNAs in treeINFERRED["CNAs"]:
        for CNA in CNAs:
            inferred_CNAs.append(CNA)
    true_positives = 0
    false_negatives = 0
    for CNA in true_CNAs:
        if CNA in inferred_CNAs:
            true_positives+=1
        else:
            false_negatives+=1
    return (true_positives,false_negatives)

def is_SNVsupported(node,tree,nodes_genotypes,setNodes):
    # A node is SNV-supported if it (or one of its descendants) contains a SNV, a CNLOH or a loss resulting in a LOH.
    global SNV_to_region
    if len(tree["SNVs"][node])>0:
        setNodes.add(node)
    for CNA in tree["CNAs"][node]:
        region,sign = CNA
        for locus in range(len(tree["SNVs"])):
            if SNV_to_region[locus]==region:
                if sign <=0 and locus in nodes_genotypes[node]: # CNA results in LOH
                    setNodes.add(node)
    # Check if descendants are SNV-supported
    for n in range(len(tree["SNVs"])):
        if tree["parents"][n]==node:
            is_SNVsupported(n,tree,nodes_genotypes,setNodes)
            if n in setNodes:
                setNodes.add(node)

def rec_genotypes(node,tree,nodes_genotypes):
    for SNV in tree["SNVs"][node]:
        nodes_genotypes[node].add(SNV)
    for n in range(len(tree["SNVs"])):
        if tree["parents"][n]==node:
            nodes_genotypes[n] = nodes_genotypes[node].copy()
            rec_genotypes(n,tree,nodes_genotypes)

def count_TP_FN_supported(treeTRUE,treeINFERRED):
    nodes_genotypes = [set() for x in treeTRUE["parents"]] # for each node, list all of the SNVs that are in the node, or above it
    rec_genotypes(0,treeTRUE,nodes_genotypes)
    nodes_SNVsupported= set()
    is_SNVsupported(0,treeTRUE,nodes_genotypes,nodes_SNVsupported)
        
    true_CNAs=[]
    for node in nodes_SNVsupported:
        for CNA in treeTRUE["CNAs"][node]:
            true_CNAs.append(CNA)
    inferred_CNAs=[]
    for CNAs in treeINFERRED["CNAs"]:
        for CNA in CNAs:
            inferred_CNAs.append(CNA)
    true_positives = 0
    false_negatives = 0
    for CNA in true_CNAs:
        if CNA in inferred_CNAs:
            true_positives+=1
        else:
            false_negatives+=1
    return (true_positives,false_negatives)

def count_TP_FN_unsupported(treeTRUE,treeINFERRED):
    """Only count CNAs which are located in a node which does not contain (in itself or its descendants) any SNV or CNA resulting in a LOH"""
    nodes_genotypes = [set() for x in treeTRUE["parents"]] # for each node, list all of the SNVs that are in the node, or above it
    rec_genotypes(0,treeTRUE,nodes_genotypes)
    nodes_SNVsupported= set()
    is_SNVsupported(0,treeTRUE,nodes_genotypes,nodes_SNVsupported)
    
        
    true_CNAs=[]
    for node in range(len(treeTRUE["parents"])):
        if not node in nodes_SNVsupported:
            for CNA in treeTRUE["CNAs"][node]:
                true_CNAs.append(CNA)
    inferred_CNAs=[]
    for CNAs in treeINFERRED["CNAs"]:
        for CNA in CNAs:
            inferred_CNAs.append(CNA)
    true_positives = 0
    false_negatives = 0
    for CNA in true_CNAs:
        if CNA in inferred_CNAs:
            true_positives+=1
        else:
            false_negatives+=1
    return (true_positives,false_negatives)
    
    


def compute_error_rates(args):
    false_positives=0
    true_negatives = 0
    true_positives=0
    false_negatives=0
    true_positives_supported=0
    false_negatives_supported=0
    true_positives_unsupported=0
    false_negatives_unsupported=0
    Bfalse_positives=0
    Btrue_negatives = 0
    Btrue_positives=0
    Bfalse_negatives=0
    Btrue_positives_supported=0
    Bfalse_negatives_supported=0
    Btrue_positives_unsupported=0
    Bfalse_negatives_unsupported=0
    for s in range(args.nseeds):
        tree_true = read_tree(args.i+str(s)+"seed_treeTRUE.gv")
        tree_COMPASS = read_tree(args.compass+str(s)+"seed_tree.gv")
        tree_BITSC2 = read_tree(args.bitsc2+str(s)+"seed_tree.gv")
        global SNV_to_region
        SNV_to_region = map_SNV_regions(args.i+str(s)+"seed_variants.csv")
        FP,TN = count_FP_TN(tree_true,tree_COMPASS)
        false_positives+=FP
        true_negatives+=TN
        TP,FN = count_TP_FN(tree_true,tree_COMPASS)
        true_positives+=TP
        false_negatives+=FN
        TP,FN = count_TP_FN_supported(tree_true,tree_COMPASS)
        true_positives_supported+=TP
        false_negatives_supported+=FN
        TP,FN = count_TP_FN_unsupported(tree_true,tree_COMPASS)
        true_positives_unsupported+=TP
        false_negatives_unsupported+=FN
        FP,TN = count_FP_TN(tree_true,tree_BITSC2)
        Bfalse_positives+=FP
        Btrue_negatives+=TN
        TP,FN = count_TP_FN(tree_true,tree_BITSC2)
        Btrue_positives+=TP
        Bfalse_negatives+=FN
        TP,FN = count_TP_FN_supported(tree_true,tree_BITSC2)
        Btrue_positives_supported+=TP
        Bfalse_negatives_supported+=FN
        TP,FN = count_TP_FN_unsupported(tree_true,tree_BITSC2)
        Btrue_positives_unsupported+=TP
        Bfalse_negatives_unsupported+=FN
    
    with open(args.o,"w") as outfile:
        outfile.write("FPR COMPASS " + str(false_positives / (false_positives + true_negatives)) + "\n")
        outfile.write("FNR COMPASS " + str(false_negatives / (false_negatives + true_positives)) + "\n")
        outfile.write("FNR supported COMPASS " + str(false_negatives_supported / (false_negatives_supported + true_positives_supported)) + "\n")
        if (false_negatives_unsupported + true_positives_unsupported) >0:
            outfile.write("FNR unsupported COMPASS " + str(false_negatives_unsupported / (false_negatives_unsupported + true_positives_unsupported)) + "\n")

        outfile.write("FPR BITSC2 " + str(Bfalse_positives / (Bfalse_positives + Btrue_negatives)) + "\n")
        outfile.write("FNR BITSC2 " + str(Bfalse_negatives / (Bfalse_negatives + Btrue_positives)) + "\n")
        outfile.write("FNR supported BITSC2 " + str(Bfalse_negatives_supported / (Bfalse_negatives_supported + Btrue_positives_supported)) + "\n")
        if (Bfalse_negatives_unsupported + Btrue_positives_unsupported):
            outfile.write("FNR unsupported BITSC2 " + str(Bfalse_negatives_unsupported / (Bfalse_negatives_unsupported + Btrue_positives_unsupported)) + "\n")

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type = str, help='Input basename')
    parser.add_argument('--compass', type = str, help='COMPASS trees basename')
    parser.add_argument('--bitsc2', type = str, help='BITSC2 trees basename')
    parser.add_argument('--nseeds', type = int, help='Number of random seeds')
    parser.add_argument('-o', type = str, help='Output file')
    args = parser.parse_args()
    compute_error_rates(args)
    
