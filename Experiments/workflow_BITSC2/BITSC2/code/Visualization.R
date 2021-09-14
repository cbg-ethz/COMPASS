

####################################################################
############ visualization function  ####################
####################################################################

Fit_visual <- function(targ_folder,X,D)
{
  cat("Visualizing samples in folder: ",targ_folder,"\n",sep='')
  ## find Rdata file names
  rdatas <- list.files(targ_folder)
  rdatas <- rdatas[grep('K.*Rdata',rdatas)]
  
  if(length(colnames(X)) > 0){
    samp_names <- colnames(X)
  } else {
    samp_names <- paste('cell',c(1:ncol(X)),sep='')
  }
  
  
  ### analyse each  Rdata file
  l <- length(rdatas)
  compare_prob <- list()
  
  for(i in 1:l)
  {
    cat("Handling file:",i,"\n")
    ### load rdata in this iteration
    cur_rdata <- paste(targ_folder,"/",rdatas[i],sep='')
    
    ### initiate pdf
    temp <- strsplit(rdatas[i],".",fixed=T)[[1]][1]
    cur_pdf <- paste(targ_folder,"/",temp,'_fit.pdf',sep='')
    
    
    pdf(cur_pdf,width = 8, height = 6)
    load(cur_rdata)
    
    # focus on chain 1    
    chain1 <- sample_organize(Trace)
    
    #likelihood trace
    p_vec <- get_samp(chain1,'likelihood')
    compare_prob[[i]] <- p_vec
    names(compare_prob)[i] <- paste('K',Params$K,sep='')
    
    # handling multiple trees
    all_tree <- out_trees(chain1)
    if(nrow(all_tree)>1){
      locs <- apply(all_tree,1,function(x){onetree_samp(chain1,x)})  # samples for each tree
      if(length(locs)>3) # only keep top 3 most frequent trees
      { temp <- order(sapply(locs,length),decreasing = T)[1:3]
      locs <- locs[temp]
      } 
    } else{locs=list(c(1:MCMC_par$Nsamp))}
    
    # for one tree, plot median estimation
    for(j in 1:length(locs))
    {
      #extract samples for one tree
      chain2 <- chain1[locs[[j]]] 
      cur_tree <- chain2[[1]]$Ttree
      
      #igraph: plot tree
      tree_plot(cur_tree)
      
      # median L,Z
      temp <- get_samp(chain2,'L')
      med_L <- apply(temp,c(1,2),median)
      p <- CN_plot(med_L,titlename = paste('median CNV','K=',Params$K)) +
        theme(text = element_text(size=12), plot.title = element_text(hjust = 0.5))
      plot(p)
      
      temp_Z <- get_samp(chain2,'Z')
      med_Z <- apply(temp_Z,c(1,2),median)
      p <- genotype_plot(med_Z,titlename = paste('median mutation,','K=',Params$K)) +
        theme(text = element_text(size=12), plot.title = element_text(hjust = 0.5))
      plot(p)
    }
    
    dev.off() 
  }
  
  #plot likelihood for each K
  pdf_comp <- paste(targ_folder,"/",'likelihood.pdf',sep='')
  pdf(pdf_comp)
  temp <- as.data.frame(compare_prob)
  temp <- gather(temp,key = 'K',value = 'likelihood',1:l)
  ylim <- quantile(temp$likelihood,c(0.02,0.98))
  p <- ggplot(data=temp)+geom_boxplot(aes(x=K,y=likelihood,fill=K)) +
    coord_cartesian(ylim = ylim)
  plot(p)
  
  dev.off()
}




