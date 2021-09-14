# COMPASS

COpy number and Mutations Phylogeny from Amplicon Single-cell Sequencing

This tool can be used to infer a tree of somatic events (mutations, copy number variants and copy-neutral loss of heterozygosity) that occurred in a tumor. It is specifically designed to be used for MissionBio's Tapestri data, where a small number of amplicons (50-300) are sequenced for thousands of single-cells.

## Compilation

A makefile is provided.

## Usage

`./COMPASS -i [sample_name] -o [output_name] --nchains 4 --chainlength 5000 --CNV 1 --sex female`

Where:
* -i is the input sample name, see below for the format of the input files
* -o is the output name. The output is a tree in graphviz format.
* --nchains indicates the number of MCMC chains to run in parallel
* --chainlength indicate the number of iterations in each MCMC
* --CNV can be set to 1 to use CNV, or 0 to only use mutations and CNLOH events
* --sex can be female (default, 2 X chromosomes) or male (1 X chromosome)

## Input
COMPASS takes as input 2 files:
* [sample_name]_variants.csv: each line corresponds to a variant. The first columns contain metadata and the remaining columns contain the counts of reference reads and alternative reads, separated by ":", for each cell.
* [sample_name]_regions.csv: each line corresponds to a region (typically, a gene). The first column is CHR_REGIONNAME, and the remaining columns contain the number of reads in this region, for each cell. This file is only required in case CNVs are used (--CNV 1).

The `data` directory contains an example synthetic input.
