
rm(list=ls())
################################################################
############# load packages ##################################
################################################################

#setwd('~/BITSC2')

require(tidyr)
require(copynumber)
require(ggplot2)
require(reshape)
require(dplyr)
require(coda)
require(gtools)
require(Rcpp)
require(RcppArmadillo)
require(shape)
require("igraph")
require(mclust)
require(vegan)
require(TailRank)


########################################################
######## load functions ################################
########################################################


source('BITSC2/code/assist_fun.R')
source('BITSC2/code/par_samp.R')
source('BITSC2/code/tree_samp_fun.R')
source('BITSC2/code/main_fun.R')
sourceCpp("BITSC2/code/params_mutipsi_likelihood.cpp")


#################################################################
########## MODEL INPUT ############################################
#################################################################

#Command line argument: input and output
args = commandArgs(trailingOnly=TRUE)
# Args: Input, output, sequencing depth, nsamples, n_clones

myseed <-  1               # set random seed
foldername <-  args[2]         # set output foldername
dir.create(foldername)  # folder where outputs are saved
D <- data.matrix(read.table(paste(args[1],"_DP.csv",sep=""),sep=","))
colnames(D) <- NULL
X <- data.matrix(read.table(paste(args[1],"_AD.csv",sep=""),sep=","))
colnames(X) <- NULL
psi <- rep(strtoi(args[3]),dim(D)[2]) #sequencing depth
# Genomic segments: group SNVs that are in the same region.
segments <- data.matrix(read.table(paste(args[1],"_segments.csv",sep=""),sep=","))
colnames(segments) <- NULL
#segments <- array(c(1:dim(D)[1],1:dim(D)[1]),dim=c(dim(D)[1],2))


##############################################
######## load parameter file ################
############################################
source('BITSC2/code/specify_pars.R')
MCMC_par$burnin <- strtoi(args[4])  # burnin sample size
MCMC_par$Nsamp <- 500   # number of samples for inference
MCMC_par$Ntune <- strtoi(args[4])  # number of samples used for adaptive parameter tuning
Nclone = c(strtoi(args[5]))
par_disp(Params, MCMC_par)


##############################################
######## sampling ##########################
############################################
source("BITSC2/code/sampler.R")


##################################################
########## model selection:   ###############
##################################################

#source("BITSC2/code/Model_select.R")
#pdf(paste(foldername,"/","selection.pdf",sep=""))
#BIC = BIC_fun(X,D,foldername)
#dev.off()

## Find best K
#best_K=-1
#best_BIC=100000000
#for (i in 1:length(BIC$K)){
#  if (BIC$BIC[i]<best_BIC){
#    best_BIC = BIC$BIC[i]
#    best_K = BIC$K[i]
#  }
#}



########################################################
########## get point estimates ###############
########################################################
source("BITSC2/code/point_estimate.R")
# # specify 
sample_Rdata <- paste("seed1_K",Nclone,".Rdata",sep="")
point_est <- get_point_estimate(foldername,sample_Rdata,Params)

origins = cbind (point_est[[1]]$Zo,point_est[[1]]$Lo)
write.table(origins,paste(args[2],"_originsBITSC2.csv",sep=""),row.names=FALSE,col.names=FALSE)
write.table(point_est[[1]]$Ttree,paste(args[2],"_treeBITSC2.csv",sep=""),row.names=FALSE,col.names=FALSE)







