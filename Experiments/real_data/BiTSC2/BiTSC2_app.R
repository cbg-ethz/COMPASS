
rm(list=ls())
################################################################
############# load packages ##################################
################################################################

setwd('./BiTSC2')

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


source('./code_function/assist_fun.R')
source('./code_function/par_samp.R')
source('./code_function/tree_samp_fun.R')
source('./code_function/main_fun.R')
sourceCpp("./code_function/params_mutipsi_likelihood.cpp")


#################################################################
########## MODEL INPUT ############################################
#################################################################

myseed <-  1               # set random seed
foldername <-  "temp_out"          # set output foldername
dir.create(foldername)  # folder where outputs are saved

scdata <- readRDS('example_data.RDS')
D <- scdata$obs.reads$D_drop # total reads, M * N matrix, where row represents locus, column represents cell
X <- scdata$obs.reads$X_drop # variant reads, M * N matrix. where row represents locus, column represents cell
#segments <- NULL
segments <- scdata$segment
psi <- rep(3,dim(D)[2]) #squencing depth


##############################################
######## load parameter file ################
############################################
source('specify_pars.R')
par_disp(Params, MCMC_par)

##############################################
######## sampling ##########################
############################################
source("./code_function/sampler.R")


##################################################
########## model selection   #####################
##################################################

source("./code_function/Model_select.R")
pdf(paste(foldername,"/","selection.pdf",sep=""))
BIC_fun(X,D,foldername)
dev.off()


##################################################
########## visualization   ######################
##################################################
cat("Visualizing sampling results: \n")
source("./code_function/Visualization.R")
Fit_visual(foldername,X,D)


########################################################
########## get point estimates ###############
########################################################
source("./code_function/point_estimate.R")
# # specify 
sample_Rdata <- "seed1_K4.Rdata"
point_est <- get_point_estimate(foldername,sample_Rdata)







