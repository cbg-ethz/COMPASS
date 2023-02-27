################################################################################
#
# this file is used to specify model parameters
#
########################################################################################


########################################################
############## MODEL PARAMETERS ############################
########################################################

set.seed(myseed)

MCMC_par <- Params <- list()

#### number of cells
Params$N <-  ncol(D)

#### number of loci
Params$M <- nrow(D)

#### maximum number of copy
Params$max_CN <- 4

#### maximum number of mutant copies
Params$max_mut <- 1

##### hyper parameter of phi: theta
Params$r=1.5

##### hyper parameter of pi
Params$a_pi <- 10000
Params$b_pi <- 1

##### hyper parameter in prior of Z: prob of gaining one mutant allele
Params$zeta <- 10^-2

##### read depth parameter
Params$psi <- psi

##### the standard deviation of the Normal distribution used to propose a new overdispersion parameter in the MH step
Params$ws_sd <- 18

##### the shape parameter of the prior Gamma distribution of the overdispersion parameter w and s
Params$ws_shape <- 11

##### the rate parameter of the prior Gamma distribution of the overdispersion parameter w and s
Params$ws_rate <- 0.1

########################################################
############## get segments ############################
########################################################

if(is.null(segments)){
  cat("Genome segment information not found. Generating segments with bin size 1\n")
  segments <- cbind(c(1:dim(D)[1]),c(1:dim(D)[1]))
}

Params$segments <- segments


########################################################
############  specify MCMC parameters ##################
########################################################

MCMC_par$burnin <- 500  # burnin sample size
MCMC_par$Nsamp <- 500   # number of samples for inference
MCMC_par$Ntune <- 500  # number of samples used for adaptive parameter tuning



MCMC_par$swap_interval <- 30 # make Matroplis Hastings move in every how many samples
MCMC_par$Nchain <- 5 # number of paralel chains
MCMC_par$delta_T <- 0.35

Temperature <- seq(1,by=MCMC_par$delta_T,length.out = MCMC_par$Nchain )  # temperatures

# ######################################################
# ############## adaptive tuning parameter ##############
# ######################################################
 
adapt_par <- list(theta_tune=6)
adapt_par <- lapply(c(1:MCMC_par$Nchain),function(x){adapt_par})


######################################################
##########  sequence of candidate Ks ##########
######################################################

# candidate subclone numbers K
# need K >= 2
Nclone <- c(3:10) 

