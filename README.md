# COMPASS

COpy number and Mutations Phylogeny from Amplicon Single-cell Sequencing

This tool can be used to infer a tree of somatic events (mutations and copy number alterations) that occurred in a tumor. It is specifically designed to be used for MissionBio's Tapestri data, where a small number of amplicons (50-300) are sequenced for thousands of single-cells.

## Quick start
```
git clone https://github.com/cbg-ethz/COMPASS.git
cd COMPASS
make
./COMPASS -i data/processed_data_AML_Morita2020/AML-59-001 -o AML-59-001 --nchains 4 --chainlength 5000 --CNA 1
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


### Use with Docker
```
docker run -t -v `pwd`:`pwd` -w `pwd` esollier/compass:v1.1 COMPASS -i data/processed_data_AML_Morita2020/AML-59-001 -o AML-59-001 --nchains 4 --chainlength 5000 --CNA 1
```

## Input
COMPASS takes as input 2 files:
* [sample_name]_variants.csv: each line corresponds to a variant. The first columns contain metadata and the remaining columns contain the counts of reference reads and alternative reads, separated by ":", for each cell.
* [sample_name]_regions.csv: each line corresponds to a region (typically, a gene). The first column is CHR_REGIONNAME, and the remaining columns contain the number of reads in this region, for each cell. This file is only required in case CNVs are used (--CNV 1).

The `data` directory contains an example synthetic input, as well as two preprocessed real AML samples. The `Experiments` directory contains scripts used to preprocess the loom files generated by the Tapestri pipeline, as well as workflows used to run simulations on synthetic data.

## Output
If [output_name] ends with .gv , COMPASS will only output the tree in graphviz format, which can then be plotted. Otherwise, COMPASS will produce as output:
* [output_name]_tree.gv: tree in graphviz format
* [output_name]_tree.json: tree in json format
* [output_name]_cellAssignments.tsv: hard assignments of cells to nodes, and whether or not the cell was inferred to be a doublet (in which case the node assignment is unreliable).
* [output_name]_cellAssignmentsProbs.tsv: posterior attachment probabilities of cells to nodes
* [output_name]_nodes_genotypes.tsv: Genotype of each SNV for each node (0: no mutation; 1: heterozygous; 2: homozygous mutated)
* [output_name]_nodes_copynumbers.tsv: Copy number of each region for each node

The `data/output` directory contains an example output.





