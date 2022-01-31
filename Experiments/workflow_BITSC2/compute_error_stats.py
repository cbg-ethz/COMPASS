import argparse
import os
import numpy as np
import pandas as pd
import copy


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



def count_FP_TN(treeTRUE,treeINFERRED):
    true_CNVs=[]
    for CNVs in treeTRUE["CNVs"]:
        for CNV in CNVs:
            true_CNVs.append(CNV)
    inferred_CNVs=[]
    for CNVs in treeINFERRED["CNVs"]:
        for CNV in CNVs:
            inferred_CNVs.append(CNV)
    false_positives = 0
    true_negatives = 0
    for CNV in inferred_CNVs:
        if not CNV in true_CNVs:
            false_positives+=1
    n_regions = 0
    for SNVs in treeTRUE["muts"]:
        for SNV in SNVs:
            n_regions = max(n_regions,SNV+1)
    for i in range(n_regions):
        is_negative=True
        for CNV in true_CNVs:
            if CNV[0]==i:
                is_negative=False
        if is_negative:
            not_inferred=True
            for CNV in inferred_CNVs:
                if CNV[0]==i:
                    not_inferred=False
            if not_inferred:
                true_negatives+=1
    return (false_positives,true_negatives)

def count_TP_FN(treeTRUE,treeINFERRED):
    true_CNVs=[]
    for CNVs in treeTRUE["CNVs"]:
        for CNV in CNVs:
            true_CNVs.append(CNV)
    inferred_CNVs=[]
    for CNVs in treeINFERRED["CNVs"]:
        for CNV in CNVs:
            inferred_CNVs.append(CNV)
    true_positives = 0
    false_negatives = 0
    for CNV in true_CNVs:
        if CNV in inferred_CNVs:
            true_positives+=1
        else:
            false_negatives+=1
    return (true_positives,false_negatives)

def is_SNVsupported(node,tree,nodes_genotypes,setNodes):
    if len(tree["muts"][node])>0:
        setNodes.add(node)
    for CNV in tree["CNVs"][node]:
        locus,sign = CNV
        if sign == -1 and locus in nodes_genotypes[node]: # CNV results in LOH
            setNodes.add(node)
    
    for n in range(len(tree["muts"])):
        if tree["parents"][n]==node:
            is_SNVsupported(n,tree,nodes_genotypes,setNodes)
            if n in setNodes:
                setNodes.add(node)

def rec_genotypes(node,tree,nodes_genotypes):
    for SNV in tree["muts"][node]:
        nodes_genotypes[node].add(SNV)
    for n in range(len(tree["muts"])):
        if tree["parents"][n]==node:
            nodes_genotypes[n] = nodes_genotypes[node].copy()
            rec_genotypes(n,tree,nodes_genotypes)

def count_TP_FN_supported(treeTRUE,treeINFERRED):
    nodes_genotypes = [set() for x in treeTRUE["parents"]] # for each node, list all of the SNVs that are in the node, or above it
    rec_genotypes(0,treeTRUE,nodes_genotypes)
    nodes_SNVsupported= set()
    is_SNVsupported(0,treeTRUE,nodes_genotypes,nodes_SNVsupported)
        
    true_CNVs=[]
    for node in nodes_SNVsupported:
        for CNV in treeTRUE["CNVs"][node]:
            true_CNVs.append(CNV)
    inferred_CNVs=[]
    for CNVs in treeINFERRED["CNVs"]:
        for CNV in CNVs:
            inferred_CNVs.append(CNV)
    true_positives = 0
    false_negatives = 0
    for CNV in true_CNVs:
        if CNV in inferred_CNVs:
            true_positives+=1
        else:
            false_negatives+=1
    return (true_positives,false_negatives)

def count_TP_FN_unsupported(treeTRUE,treeINFERRED):
    """Only count CNVs which are located in a node which does not contain (in itself or its descendants) any SNV or CNV resulting in a LOH"""
    nodes_genotypes = [set() for x in treeTRUE["parents"]] # for each node, list all of the SNVs that are in the node, or above it
    rec_genotypes(0,treeTRUE,nodes_genotypes)
    nodes_SNVsupported= set()
    is_SNVsupported(0,treeTRUE,nodes_genotypes,nodes_SNVsupported)
    
        
    true_CNVs=[]
    for node in range(len(treeTRUE["parents"])):
        if not node in nodes_SNVsupported:
            for CNV in treeTRUE["CNVs"][node]:
                true_CNVs.append(CNV)
    inferred_CNVs=[]
    for CNVs in treeINFERRED["CNVs"]:
        for CNV in CNVs:
            inferred_CNVs.append(CNV)
    true_positives = 0
    false_negatives = 0
    for CNV in true_CNVs:
        if CNV in inferred_CNVs:
            true_positives+=1
        else:
            false_negatives+=1
    return (true_positives,false_negatives)
    
    


def compute_error_rates(basename,nseeds,output):
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
    for s in range(nseeds):
        try:
            tree_true = read_tree_gv(basename+"_"+str(s)+"_treeTRUE.gv")
            tree_COMPASS = read_tree_gv(basename+"_"+str(s)+"_treeCOMPASS.gv")
            tree_BITSC2 = read_tree_BITSC2(basename+"_"+str(s))
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
        except:
            print("failed "+str(s))
    
    with open(output,"w") as outfile:
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
    parser.add_argument('--nseeds', type = int, help='Number of random seeds')
    parser.add_argument('-o', type = str, help='Output file')
    args = parser.parse_args()
    compute_error_rates(args.i,args.nseeds,args.o)
    