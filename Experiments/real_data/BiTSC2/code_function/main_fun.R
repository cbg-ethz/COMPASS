##################################################
##### parameter random initiation ###############
##################################################

init_gen_II <- function(Params, Heu_tree_cluster=NULL)
{ # psi: read depth
  
  out <- list()
  
  
  #initial tree
  if(is.null(Heu_tree_cluster)){
    out$Ttree <- c(0, 1, rep(1,Params$K-2))
    if(Params$K>2){
      for(z in 3:Params$K){out$Ttree[z] <- sample(1:(z-1),1)}
    }
  } else {
    out$Ttree <- Heu_tree_cluster$tree
  }
  
  # initial pi
  out$pi <- 0.5
  
  # initial rho
  out$rho <- runif(1,min = 0, max = 1)
  
  # initial mu
  out$mu <- runif(1,min = 0, max = 1)
  
  # initial w
  out$w <- rgamma(1,shape = Params$ws_shape, rate  = Params$ws_rate)
  
  # initial s
  out$s <- rgamma(1,shape = Params$ws_shape, rate  = Params$ws_rate)
  
  out$g <- rep(0,Params$M)
  # initial L
  tempL <- matrix(1, Params$M, 2)
  tempZ <- matrix(1, Params$M, Params$K)
  segs <- Params$segments
  for(i in 1:nrow(segs))
  {
    cur_seg <- seq(segs[i,1], segs[i,2])
    temp <- propose_Lo_cpp(out$Ttree, tempZ[cur_seg,,drop=F], Params)
    for(j in cur_seg){
      tempL[j,] <- temp
    }
  }
  out$Lo <- tempL
  out$L <- Lo_to_L_cpp(out$Lo, out$Ttree)
  
  # initial Z
  tempZ <- matrix(0, Params$M, 2)
  for(j in 1:Params$M)
  {
    tempZ[j,] <- propose_Zo_cpp(out$Ttree,out$L[j,],Params)
  }
  out$Zo <- tempZ
  out$Z <- Zo_to_Z_cpp(out$Zo,out$Lo,out$g,out$Ttree)
  
  # initial theta
  out$theta <- Heu_tree_cluster$phi*10
  
  # initial phi
  out$phi <- Heu_tree_cluster$phi
  
  # initial C
  out$C <- Heu_tree_cluster$cluster
  
  out
}



#######################################################
########## heuristic init gen ####################
#######################################################
heuristic_tree_cluster_gen <- function(X,D,K)
{
  
  P <- t(X/D) # observed VAF
  P[which(t(D)==0)] <- 0
  P_cluster <- Mclust(P,G=K)
  centers <- t(P_cluster$parameters$mean)
  vaf_mean <- apply(centers,1,mean)
  ord <- order(vaf_mean,decreasing = F)
  centers <- centers[ord,]
  
  # initialize cluster
  C <- P_cluster$classification
  C <- sapply(C,function(x){which(ord==x)})
  
  # initialize phi
  phi <- P_cluster$parameters$pro[ord]
  
  # initialize tree
  dis<- dist(centers,p=2)
  dis <- as.matrix(dis)
  spanningtree <- spantree(dis)
  rootid <- which.min(dis[2:nrow(dis),1])
  idspantree<-spanningtree$kid
  tree0 <- c(0,idspantree)
  
  ##### sort new tree if violate order constraint
  temp <- tree_org(tree0)
  new_tree <- temp$tree
  new_phi <- phi
  new_C <- C
  
  if(!all(temp$ord==c(1:K))){
    new_phi <- new_phi[temp$ord]
    new_C <- sapply(C,function(x){which(temp$ord==x)})
  }

  tree_cluster <- list(tree=new_tree, cluster=new_C, phi=new_phi)
  return(tree_cluster)
}



#######################################################
########## acquire one MCMC sample ####################
#######################################################

MCMC_onesamp_II <- function(X, D, Params, init, adap, temper=1)
{ # init: last round MCMC sample
  # temper: temperature
  # adap: adaptive tuning parameters
  
  out <- init
  
  ########## update rho #######
  
  out$rho <- samp_rho(out,Params,D,temper)
  
  ########## update pi #######
  
  out$pi <- samp_pi(out,Params,temper)
  
  ########## update w #######
  
  out$w <- samp_w(out,Params,X,D,temper)
  
  ########## update s #######
  
  out$s <- samp_s(out,Params,D,temper)
  
  ########## update mu #######
  
  out$mu <- samp_mu(out,Params,X,D,temper)
  
  ######### update L matrix ##########
  
  out$Lo <- samp_Lo_all(out,X,D,Params,temper)
  
  out$L <- Lo_to_L_cpp(out$Lo, out$Ttree)
  
  #### update Z matrix #########
  
  out$Zo <- samp_Zo_all(out, X, D, Params, temper)
  
  out$g <- estimate_g(out, X, D, Params)
  
  out$Z <- Zo_to_Z_cpp(out$Zo, out$Lo, out$g, out$Ttree)
  
  ##### update tree structure ########
  
  if(Params$K>2)
  {
    Lo <- out$Lo
    Zo <- out$Zo
    
    u <- runif(1)
    if(u<0.1){
      slice_out <- slice_samp_tree(X,D,out,Lo,Zo,Params,temper=temper)
      if(!is.null(slice_out))
      {
        out$Ttree <- slice_out$Ttree
        out$L <- slice_out$L
        out$Z <- slice_out$Z
      }
    } else {
      out <- MH_samp_tree(out$Ttree,X,D,out,Params,Lo,Zo,temper)
    }
  }
  
  ########## update phi (or Theta)###
  
  out$theta <- samp_theta_all(out,Params,adap$theta_tune,temper)
  out$phi <- out$theta/sum(out$theta) # update phi with new theta values
  
  
  ########## update C ########
  
  out$C <- samp_C_all(out,Params,X,D,temper)
  
  
  ###### keep track of likelihood ########
  #p=prob_XD(X,D,out$phi,out$L,out$Z,Params$psi)
  LL <- out$L[,out$C]
  ZZ <- out$Z[,out$C]
  p <- sum(log_prob_X(X,D,LL,ZZ,out$w,out$mu)+log_prob_D(D,LL,Params$psi,out$rho,out$s,out$mu))
  out$likelihood <- p
  
  out
}



















