################################################################
# this file contains sampling procedure procedure of model
##################################################################



for( Nc in Nclone ){
  if(Nc<2) stop("K need to be >= 2")
  
  ##### number of subclones
  Params$K <- Nc
  
  
  ######################################################
  ######## initial value for parameters  ###############
  ######################################################
  
  New_par <- vector('list',MCMC_par$Nchain)
  
  for(i in 1:MCMC_par$Nchain){
    if(Params$K==2){
      tree0 <- c(0,1)
    } else {#tree0 <- heuristic_tree_gen(X,D,Params$K)
      HTC <- heuristic_tree_cluster_gen(X,D,Params$K)
      }
    #New_par[[i]] <- init_gen(Params,tree0)
    New_par[[i]] <- init_gen_II(Params,HTC)
  }

  #################################################################
  ################### MCMC sampling for each K #####################
  ################################################################
  
  start_t <- Sys.time()
  
  Nrep <- MCMC_par$burnin+ MCMC_par$Nsamp + MCMC_par$Ntune # number of total samples
  
  tune_samp <- vector('list',MCMC_par$Nchain)
  Trace <- vector('list',MCMC_par$Nsamp)
  beta_log_like <- matrix(0,MCMC_par$Nchain ,MCMC_par$Nsamp)  # keep beta chain likelihood for WBIC
  # Trace <- vector('list',Nrep)
  # beta_log_like <- matrix(0,MCMC_par$Nchain, Nrep)
  
  for(i in 1:MCMC_par$Nchain) {
    tune_samp[[i]] <- vector('list',MCMC_par$Ntune)
  }
  
  for(i in 1:Nrep){
    
    ######## perform MCMC on each chain ##########
    for (z in 1:MCMC_par$Nchain){
      temp <- MCMC_onesamp_II(X,D,Params,New_par[[z]],adapt_par[[z]],temper=Temperature[z])
      New_par[[z]] <- temp
    }
    
    
    #####  tuning phase  ########
    
    if(i <= MCMC_par$Ntune){
      for(z in 1:MCMC_par$Nchain){
        tune_samp[[z]][[i]] <- New_par[[z]]
        
        if( i %% 50==0) # check mixing at certain intevals
        {
          sub=c((i-49):i)
          
          theta_100 <- get_samp(tune_samp[[z]],'theta',sub)  # previous 100 theta samples
          
          # theta
          temp <- adaptive_tune(theta_100,low=0.4,high=0.65)  # determine how to change tuning parameter
          if(temp$move==1)
          {adapt_par[[z]]$theta_tune <- adapt_par[[z]]$theta_tune+1} else if(temp$move==-1)
          {adapt_par[[z]]$theta_tune <- max(adapt_par[[z]]$theta_tune-1,1)}
          
        }
      }
    }
    
    
    ############### switch between chains ###########################
    
    if(MCMC_par$Nchain>1){
      if(i%%MCMC_par$swap_interval == 0)
      {
        a <- sample(c(1:(MCMC_par$Nchain-1)),1) # randomly select a chain
        b <- a+1
        #log  of numerater
        logp_num <- New_par[[a]]$likelihood/Temperature[b] + New_par[[b]]$likelihood/Temperature[a]  
        logp_num <-  logp_num + log_prior_all(New_par[[a]],Params,Temperature[b]) +
          log_prior_all(New_par[[b]],Params,Temperature[a])
        #log of denominator
        logp_denom <-  New_par[[b]]$likelihood/Temperature[b] + New_par[[a]]$likelihood/Temperature[a]   
        logp_denom <-  logp_denom + log_prior_all(New_par[[a]],Params,Temperature[a]) +
          log_prior_all(New_par[[b]],Params,Temperature[b])
        acc_prob <- min(1,exp(logp_num-logp_denom))  # probability of accepting swap
        
        if(runif(1)<acc_prob)
        {
          temp <- New_par[[a]]
          New_par[[a]] <- New_par[[b]]
          New_par[[b]] <- temp
        }
        
      }
    }
    
    
    ######## after burnin stage: save new samples into Trace ###############
    #### only keep chain 1
    ct <- i- MCMC_par$burnin - MCMC_par$Ntune  # count number of saved samples
    if(ct > 0) {
      Trace[[ct]] <- New_par[[1]]
      for(z in 1:length(New_par)){
        beta_log_like[z,ct] <-  New_par[[z]]$likelihood # calculate Ln for all chains
      }
    }

    
    ########## display progress ###############
    
    if( i %% 100 == 0){
      ave_speed <- difftime(Sys.time(),start_t,units = "mins")/i
      cat("Subclone number K =",Nc,"accomplished",100*i/(MCMC_par$Nsamp+MCMC_par$burnin+MCMC_par$Ntune),
          "%; time remaining: ",(MCMC_par$Nsamp+MCMC_par$burnin+MCMC_par$Ntune-i)*ave_speed,"mins \n")
    }
  }
  
  
  cat('Subclone number:',Nc,'completed. Time consumed:',difftime(Sys.time(),start_t,units = "mins"),'mins \n \n')
  
  
  ##################################################################
  ###### save data ############################################
  ##################################################################
  
  cur_file=paste(foldername,'/seed',myseed,'_K',Nc,'.Rdata',sep='')
  save(Params,MCMC_par,Trace,myseed,beta_log_like,file=cur_file)
}
