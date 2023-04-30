# COMPASS

COpy number and Mutations Phylogeny from Amplicon Single-cell Sequencing

This tool can be used to infer a tree of somatic events (mutations and copy number alterations) that occurred in a tumor. It is specifically designed to be used for MissionBio's Tapestri data, where a small number of amplicons (50-300) are sequenced for thousands of single-cells.

## Quick start
```
git clone https://github.com/cbg-ethz/COMPASS.git
cd COMPASS
make
./COMPASS -i data/preprocessed_data_AML_Morita2020/AML-59-001 -o AML-59-001 --nchains 4 --chainlength 5000 --CNA 1
dot -Tpng -o AML-59-001_tree.png AML-59-001_tree.gv
```

Graphviz is required in order to plot the tree, which can be installed on Ubuntu by running `sudo apt-get install graphviz `

## Usage

`./COMPASS -i [sample_name] -o [output_name] --nchains 4 --chainlength 5000 --CNA 1 --sex female`

Where:
* -i is the input sample name, see below for the format of the input files
* -o is the output name. The output is a tree in graphviz format.
* --nchains indicates the number of MCMC chains to run in parallel
* --chainlength indicate the number of iterations in each MCMC
* --CNA can be set to 1 to use CNA, or 0 to only use SNVs
* --sex can be female (default, 2 X chromosomes) or male (1 X chromosome)
Additional parameters can be changed if needed, although their default values should work for most cases:
* -d (default: 1): if 1, COMPASS will use the model with doublets, and if 0, COMPASS will use the model without doubets (faster)
* --doubletrate (default: 0.08) determines the doublet rate, in case -d is set to 1.
* --dropoutrate (default: 0.05): prior mean of the allelic dropout rates. The dropout rates will be estimated for each SNV with a beta binomial distribution.
* --dropoutrate_concentration (default: 100): prior concentration parameter for the beta binomial distribution for dropout rates. Higher values will result in the estimated dropout rates to be closer to the prior mean.
* --seqerror (default: 0.02): sequencing error rate
* --nodecost (default: 1): cost of adding a node to the tree
* --cnacost (default: 85): cost of adding a CNA event to the tree
* --lohcost (default: 85) cost of adding a LOH event to the tree

In targeted sequencing, different regions have different coverages, depending on the number of amplicons targeting each region and the efficiency of the primers. By default, COMPASS will use the cells attached to the root in order to estimate the proportion of reads falling on each region in the absence of CNAs. Optionally, it is possible to provide the weights of each region with the arguments --regionweights. An example csv file is provided in `data/preprocessed_data_Morita2020/region_weights_50amplicons.csv` and a script to generate such a csv file is provided at `Experiments/preprocessing/estimate_region_weights.py`.


### Use with Docker
```
docker run -t -v `pwd`:`pwd` -w `pwd` esollier/compass:v1.1 COMPASS -i data/preprocessed_data_AML_Morita2020/AML-59-001 -o AML-59-001 --nchains 4 --chainlength 5000 --CNA 1
```

## Input
COMPASS takes as input 2 files:
* [sample_name]_variants.csv: each line corresponds to a variant. The first columns contain metadata and the remaining columns contain the counts of reference reads and alternative reads, separated by ":", for each cell.
* [sample_name]_regions.csv: each line corresponds to a region (typically, a gene). The first column is CHR_REGIONNAME, and the remaining columns contain the number of reads in this region, for each cell. This file is only required in case CNVs are used (--CNV 1).

The `data` directory contains preprocessed datasets. The `Experiments/preprocessing` directory contains scripts used to preprocess the loom files generated by the Tapestri pipeline, as well as workflows used to run simulations on synthetic data.

COMPASS contains several hardcoded parameters which are defined at the end of `input.cpp`. We expect these parameters to be robust and that users will not have to adapt them, but if necessary, it is possible to change them and recompile COMPASS. For example, increasing or decreasing the CNA_cost parameter will result in trees with less/more CNAs, or increasing/decreasing the node_cost parameter will result in trees with fewer/more nodes.

## Output
If [output_name] ends with .gv , COMPASS will only output the tree in graphviz format, which can then be plotted. Otherwise, COMPASS will produce as output:
* [output_name]_tree.gv: tree in graphviz format
* [output_name]_tree.json: tree in json format
* [output_name]_cellAssignments.tsv: hard assignments of cells to nodes, and whether or not the cell was inferred to be a doublet (in which case the node assignment is unreliable).
* [output_name]_cellAssignmentsProbs.tsv: posterior attachment probabilities of cells to nodes
* [output_name]_nodes_genotypes.tsv: Genotype of each SNV for each node (0: no mutation; 1: heterozygous; 2: homozygous mutated)
* [output_name]_nodes_copynumbers.tsv: Copy number of each region for each node

The `data/output_example` directory contains an example output.





