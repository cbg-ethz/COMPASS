####################################################################
############ get point estimation  ####################
####################################################################

get_point_estimate <- function(foldername,rdata,Params)
{
  ### load rdata in this iteration
  cur_rdata <- paste(foldername,rdata,sep='/')
  load(cur_rdata)
  # reorganize samples
  chain1 <- sample_organize(Trace)
  # handling multiple trees
  all_tree <- out_trees(chain1)
  if(nrow(all_tree)>1){
    locs <- apply(all_tree,1,function(x){onetree_samp(chain1,x)})  # samples for each tree
    if(length(locs)>3) # only keep top 3 most frequent trees
    { temp <- order(sapply(locs,length),decreasing = T)[1:3]
    locs <- locs[temp]
    } 
  } else{locs <- list(c(1:MCMC_par$Nsamp))}
    
  # for one tree, plot median estimation
  out <- list()
  for(j in 1:length(locs))
  {
    point_est <- list()
    #extract samples for one tree
    chain2 <- chain1[locs[[j]]] 
    cur_tree <- chain2[[1]]$Ttree
    point_est$Ttree <- cur_tree

    # mean phi
    temp <- get_samp(chain2,'phi')
    med_phi <- apply(temp,2,mean)
    point_est$phi <- med_phi
    
    # mean rho
    temp <- get_samp(chain2,'rho')
    med_rho <- mean(temp)
    point_est$rho <- med_rho
    
    # mean mu
    temp <- get_samp(chain2,'mu')
    med_mu <- mean(temp)
    point_est$mu <- med_mu

    # mean w
    temp <- get_samp(chain2,'w')
    med_w <- mean(temp)
    point_est$w <- med_w

    # mean s
    temp <- get_samp(chain2,'s')
    med_s <- mean(temp)
    point_est$s <- med_s
    
    # mean C
    temp <- get_samp(chain2,'C')
    med_C <- apply(temp,2,Mode)
    point_est$C <- med_C
    
    # median L
    temp <- get_samp(chain2,'L')
    med_L <- apply(temp,c(1,2),median)
    point_est$L <- med_L
    
    # median Z
    temp_Z <- get_samp(chain2,'Z')
    med_Z <- apply(temp_Z,c(1,2),median)
    point_est$Z <- med_Z
    
    
    # estimation of g
    Zo <- Z_to_Zo_cpp(point_est$Z)
    Lo <- L_to_Lo_cpp(point_est$L)
    point_est$Zo <- Zo
    point_est$Lo <- Lo
    point_est$g <- estimate_g(point_est,X,D,Params)
    
    out[[j]] <- point_est
  }
  return(out)
}


Mode <- function(x) {
  ux <-  unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

