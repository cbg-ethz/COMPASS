
rm(list=ls())
################################################################
############# load packages ##################################
################################################################

#setwd('./BiTSC2')

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
require(igraph)
require(mclust)
require(vegan)
require(TailRank)


########################################################
######## load functions ################################
########################################################


source('BiTSC2/code_function/assist_fun.R')
source('BiTSC2/code_function/par_samp.R')
source('BiTSC2/code_function/tree_samp_fun.R')
source('BiTSC2/code_function/main_fun.R')
sourceCpp("BiTSC2/code_function/params_mutipsi_likelihood.cpp")


#################################################################
########## MODEL INPUT ############################################
#################################################################

# Args: Input, output, nsamples, sequencing depth , n_clones
args = commandArgs(trailingOnly=TRUE)

myseed <-  1               # set random seed
foldername <-  args[2]          # set output foldername
dir.create(foldername)  # folder where outputs are saved

D <- data.matrix(read.table(paste(args[1],"_DP.csv",sep=""),sep=","))
colnames(D) <- NULL
X <- data.matrix(read.table(paste(args[1],"_AD.csv",sep=""),sep=","))
colnames(X) <- NULL



# Genomic segments: group SNVs that are in the same region.
segments <- data.matrix(read.table(paste(args[1],"_segments.csv",sep=""),sep=","))
colnames(segments) <- NULL


##############################################
######## load parameter file ################
############################################
if (length(args)<=3){
  psi <- rep(20,dim(D)[2])
  for (j in 1:dim(D)[2]){
	  psi[j] = mean(D[,j])
  }
}else{
  psi <- rep(strtoi(args[4]),dim(D)[2]) #sequencing depth
}
source('BiTSC2/specify_pars.R')
MCMC_par$burnin <- strtoi(args[3])  # burnin sample size
MCMC_par$Nsamp <- 500   # number of samples for inference
MCMC_par$Ntune <- strtoi(args[3])  # number of samples used for adaptive parameter tuning

if (length(args)<=3){
  Nclone = c(3,4,5,6,7) #3,4,5,6
}else{
  Nclone = c(strtoi(args[5])) 
}


par_disp(Params, MCMC_par)

##############################################
######## sampling ##########################
############################################
source("BiTSC2/code_function/sampler.R")


##################################################
########## model selection   #####################
##################################################

source("BiTSC2/code_function/Model_select.R")
pdf(paste(foldername,"/","selection.pdf",sep=""))
BIC <- BIC_fun(X,D,foldername)
dev.off()

## Find best K
best_K=-1
best_BIC=100000000
for (i in 1:length(BIC$K)){
  if (BIC$BIC[i]<best_BIC){
    best_BIC = BIC$BIC[i]
    best_K = BIC$K[i]
  }
}


##################################################
########## visualization   ######################
##################################################
cat("Visualizing sampling results: \n")
source("BiTSC2/code_function/Visualization.R")
Fit_visual(foldername,X,D)


########################################################
########## get point estimates ###############
########################################################
source("BiTSC2/code_function/point_estimate.R")
# # specify 
sample_Rdata <- paste("seed1_K",best_K,".Rdata",sep="")
point_est <- get_point_estimate(foldername,sample_Rdata)

origins = cbind (point_est[[1]]$Zo,point_est[[1]]$Lo,point_est[[1]]$g)
write.table(origins,paste(args[2],"_originsBITSC2.csv",sep=""),row.names=FALSE,col.names=FALSE)
write.table(point_est[[1]]$Ttree,paste(args[2],"_treeBITSC2.csv",sep=""),row.names=FALSE,col.names=FALSE)

attachments = point_est[[1]]$C
write.table(attachments,paste(args[2],"_attachmentsBITSC2.csv",sep=""),row.names=FALSE,col.names=FALSE)






