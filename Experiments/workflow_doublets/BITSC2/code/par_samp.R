
#####################################################
########## sampling of pi #########################
#######################################################

samp_pi <- function(samp, Params, temper)
{
  L <- samp$L
  temp <- apply(L,1,function(x){all(x==2)})
  n <- sum(temp)
  m <- Params$M-n
  rbeta(1,(n+Params$a_pi+temper-1)/temper, (m+Params$b_pi+temper-1)/temper)
}

#####################################################
########## sampling of base dropout rate rho #########################
#######################################################

samp_rho <- function(samp, Params, D, temper)
{
  rho0 <- samp$rho
  rho <- rho0
  LL <- samp$L[,samp$C]
  s <- samp$s
  mu <- samp$mu
  psi <- Params$psi
  
  rho1 <- runif(1,min=0,max=1)
  p0 <- sum(log_prob_D(D, LL, psi, rho0, s, mu))/temper
  p1 <- sum(log_prob_D(D, LL, psi, rho1, s, mu))/temper
  alpha <- exp(p1-p0)
  if(runif(1,min=0,max=1)<min(1,alpha)){rho <- rho1}
  rho
}

########################################################################
########### sampling of allelic dropout rate mu #########################
########################################################################

samp_mu <- function(samp, Params, X, D, temper)
{
  LL <- samp$L[,samp$C]
  ZZ <- samp$Z[,samp$C]
  rho <- samp$rho
  w <- samp$w
  s <- samp$s
  # propose new value
  mu0 <- samp$mu
  mu <- mu0
  mu1 <- runif(1,min=0,max=1)
  
  p0 <- sum(log_prob_X(X,D,LL,ZZ,w,mu0)+log_prob_D(D, LL, psi, rho, s, mu0))/temper
  p1 <- sum(log_prob_X(X,D,LL,ZZ,w,mu1)+log_prob_D(D, LL, psi, rho, s, mu1))/temper
  
  alpha <- exp(p1-p0)
  if(runif(1,min=0,max=1)<min(1,alpha)){mu <- mu1}
  mu
}

########################################################################
############# sampling of w ############################################
########################################################################

samp_w <- function(samp, Params, X, D, temper){
  LL <- samp$L[,samp$C]
  ZZ <- samp$Z[,samp$C]
  mu <- samp$mu
  # propose new value
  w0 <- samp$w
  w1 <- rnorm(1L, mean=w0, sd=Params$ws_sd)
  # if allowed
  if(w1 > 0){
    # compute priors
    lpriors <- dgamma(x = c(w0, w1),shape = Params$ws_shape, rate  = Params$ws_rate, log = TRUE)/temper
    
    # compute loglikelihoods
    p0 <- sum(log_prob_X(X,D,LL,ZZ,w0,mu))/temper
    p1 <- sum(log_prob_X(X,D,LL,ZZ,w1,mu))/temper
    
    # evaluate MH ratio
    if(runif(n=1L, min=0, max=1) <= exp(p1 + lpriors[2] - p0 - lpriors[1])){
      return(w1)
    }else{
      return(w0)
    }
  }else{
    return(w0)
  }
}


########################################################################
############# sampling of s ############################################
########################################################################

samp_s <- function(samp, Params, D, temper){
  LL <- samp$L[,samp$C]
  psi <- Params$psi
  rho <- samp$rho
  mu <- samp$mu
  
  # propose new value
  s0 <- samp$s
  s1 <- rnorm(1L, mean=s0, sd=Params$ws_sd)
  # if allowed
  if(s1 > 0){
    # compute priors
    lpriors <- dgamma(x = c(s0, s1),shape = Params$ws_shape, rate  = Params$ws_rate, log = TRUE)/temper
    
    # compute loglikelihoods
    p0 <- sum(log_prob_D(D, LL, psi, rho, s0, mu))/temper
    p1 <- sum(log_prob_D(D, LL, psi, rho, s1, mu))/temper
    
    # evaluate MH ratio
    if(runif(n=1L, min=0, max=1) <= exp(p1 + lpriors[2] - p0 - lpriors[1])){
      return(s1)
    }else{
      return(s0)
    }
  }else{
    return(s0)
  }
}



##################################################################
########## tune theta parameter ###############################
##################################################################
adaptive_tune <- function(S,low=0.4,high=0.6)
{
  mcmcS <-  mcmc(S)
  rej_rate <- rejectionRate(mcmcS)
  rej_rate <- mean(rej_rate)
  
  if(rej_rate>high)
  {out <- 1} else if (rej_rate<low)
  {out <- -1} else {out <- 0}
  out2 <- list(move=out,rate=rej_rate)
  out2
}


