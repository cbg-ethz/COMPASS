# Estimate region weights using samples that do not have copy number alterations (if a CNA is present in less than 50% of the samples it is still OK).
# Usage: python estimate_region_weights.py output.csv sample1_regions.csv sample2_regions.csv sample3_regions.csv
# where the sampleX_regions.csv files contain the read counts for each region and each cell (same as COMPASS input)


import sys
import numpy as np
import pandas as pd

regions_proportions={}
output_file = sys.argv[1]
for infile in sys.argv[2:]:
    df = pd.read_csv(infile,sep=",",index_col=0,header=None)
    df = df / np.sum(df)
    for i in df.index:
        region = i.split("_")[-1]
        if not region in regions_proportions:
            regions_proportions[region] = []
        for j in df.columns:
            regions_proportions[region].append(df.loc[i,j])

regions_weights = {}
for region in regions_proportions:
    regions_weights[region] = np.median(regions_proportions[region])


with open(output_file,"w") as out:
    for region in regions_weights:
        tmp = out.write(region+","+str(regions_weights[region])+"\n")
