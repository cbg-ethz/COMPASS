#### assisting functions ############

### plot genotype matrix #####
genotype_plot <- function(Z,xname='subclone',yname='mutation',titlename=NULL)
{
 rownames(Z) <- colnames(Z) <- c()
 Z <- as.matrix(Z)

 cc <- colorRampPalette(c('#FFF6F1','#741521'))(length(unique(as.numeric(Z))))
 dp=melt(Z)
 colnames(dp) <- c("loci","subclone","Value")
 p <- ggplot(dp, aes(subclone,loci)) + geom_tile(aes(fill = as.factor(Value)))+
   theme(axis.text.x=element_text(angle = 90))+ 
   guides(fill = guide_legend(title = "State")) +scale_fill_manual(values=cc)+
   theme_bw()+xlab(xname) +ylab(yname) +
   theme(panel.grid =element_blank()) + 
   theme(axis.text.x = element_text(size = 20, face = "bold", vjust = 0.5, hjust = 0.5, colour = 'black'),axis.text.y = element_text(size = 15, face = "bold", vjust = 0.5, hjust = 0.5,colour = 'black'))+
   theme(axis.title = element_text(size = 18, face = "bold", vjust = 0.5, hjust = 0.5))+
   theme(legend.text = element_text(size = 15,face = 'bold'),legend.title = element_text(size = 15,face = 'bold'))
 p
}

### plot copy number matrix ######
CN_plot <- function(L,xname='subclone',yname='mutation',titlename=NULL)
{
  rownames(L) <- colnames(L) <- c()
  L <- as.matrix(L)
  dp=melt(L) 
  colnames(dp) <- c("loci","subclone","Value")
  cc <- colorRampPalette(c('#E2EBF3','#4D85BD','#B2473E'))(length(unique(as.numeric(L))))
  p <- ggplot(dp, aes(subclone,loci)) + geom_tile(aes(fill = as.factor(Value)))+
    theme(axis.text.x=element_text(angle = 90))+ 
    guides(fill = guide_legend(title = "Copy\nNumber")) + scale_fill_manual(values=cc)+
    theme_bw()+xlab(xname) +ylab(yname) +
    theme(panel.grid =element_blank()) + 
    theme(axis.text.x = element_text(size = 20, face = "bold", vjust = 0.5, hjust = 0.5, colour = 'black'),axis.text.y = element_text(size = 20, face = "bold", vjust = 0.5, hjust = 0.5,colour = 'black'))+
    theme(axis.title = element_text(size = 18, face = "bold", vjust = 0.5, hjust = 0.5))+
    theme(legend.text = element_text(size = 15,face = 'bold'),legend.title = element_text(size = 15,face = 'bold'))
  p
}



################################################
################ tree plot ####################
################################################

tree_plot <- function(tree)
{
  temp <- tree
  temp[1] <- 1
  tree_df <- data.frame(parent=temp,id=c(1:length(temp)))
  tree_df <- tree_df[-1,]
  tree_graph <- graph.data.frame(tree_df)
  col <- intpalette(c('magenta','red','limegreen','skyblue','pink','brown','gold','blue','cyan'),length(tree))
  plot(tree_graph,layout=layout_as_tree,vertex.color=col,edge.color='black',
       vertex.label.color='black',edge.lty=c("solid"),vertex.frame.color='white',
       edge.width=4)
}



##################################################################
############## get samples for given parameters ##################
##################################################################

get_samp <- function(MCMCout,var_name,subset=NULL)  
{
  #subset: a subset of samples
  
  if(is.null(subset)){subset <- c(1:length(MCMCout))}
  MCMCout <- MCMCout[subset]
  
  Ns <- length(MCMCout)  # number of samples
  temp <- MCMCout[[1]][[var_name]]
  if(is.vector(temp)){temp <- matrix(temp,nrow = 1)}
  if(is.matrix(temp)){
    if(ncol(temp)==1){
      temp <- t(temp)
      }
    }
  
  d1 <- nrow(temp)  
  d2 <- ncol(temp)
  
  if(is.null(d1))#no row dimension,then a scalar
  {
    out <- c(1:Ns)
    for(z in 1:Ns)
    {
      temp <- MCMCout[[z]][[var_name]]
      out[z] <- temp
    }
  } else if(d1==1)  # if parameter is vector, get samples as a matrix
  {
    out <- matrix(0,Ns,d2)
    for(z in 1:Ns)
    {
      temp <- MCMCout[[z]][[var_name]]
      out[z,] <- temp
    }
  } else {  # else get samples as an array
    out <- array(0,dim=c(d1,d2,Ns))
    for(z in 1:Ns)
    {
      temp <- MCMCout[[z]][[var_name]]
      out[,,z] <- temp
    }
  }
  out
}

################################################################
################ evaluate one posteriror sample ################
################################################################

clone_sort <- function(Z1,Z0)  # pick permutation with minimum error on Z
{
  #Z1: sample mutation matrix
  #Z0: true mutation matrix
  a <- ncol(Z0)
  M <- nrow(Z0)
  ords <- permutations(a-1,a-1,v=2:a) # all possible permutations excluding normal subclone
  ords <- cbind(1,ords)
  
  l <- nrow(ords) # number of possible orders
  
  loc <- 0
  v <- Inf
  for(i in 1:l){
    Znew <- Z1[,ords[i,]]
    temp <- sum(abs(Znew-Z0)) # measure of this permutation
    if(temp < v){
      v <- temp
      loc <- i
    }
  }
  
  v <- v/(a*M)
  
  out <- list(measure=v,ord=ords[loc,])
  out
}


############################################################################
################ get samples with different tree structures ################
############################################################################
# get all tree structures from posterior samples
out_trees <- function(MCMCout)
{
  l <- length(MCMCout)
  all_tree <- get_samp(MCMCout,'Ttree')
  
  all_tree <- all_tree[!duplicated(all_tree),,drop=F]
  
  out <- all_tree
  out
}

## get sample numbers for one tree structure
onetree_samp <- function(MCMCout,tree)
{
  l <- length(MCMCout)
  
  a <- c(1:l)
  for(i in 1:l){
    temp <- MCMCout[[i]]$Ttree
    if(all(temp==tree)){a[i] <- 1} else{a[i] <- 0}
  }
  
  out <- which(a==1)
  out
}


### organize tree by new order
get_newtree <- function(t0, new_ord){
  stopifnot(length(new_ord)==length(t0))
  l <- length(new_ord)
  new_t <- c(0:(l-1))
  if (l>2){
    for (i in 2:l){
      new_t[i] <- which(new_ord==t0[new_ord[i]])
    }
  }
  return(new_t)
}


### reduce duplicated trees
sample_organize <- function(trace) {
  # organize 
  l <- length(trace)
  stopifnot(l>=2)
  # organize samples
  Z_list <- list()
  Z_list[[1]] <- trace[[1]]$Z
  K <- ncol(trace[[1]]$Z)
  
  # re-organize each sample
  for(i in 2:length(trace)){
    cur_Z <- trace[[i]]$Z
    err_vec <- c(1:length(Z_list))
    ord_list <- list()
    for(j in 1:length(Z_list)){
      c_sort <- clone_sort(cur_Z,Z_list[[j]])
      err_vec[j] <- c_sort$measure
      ord_list[[j]] <- c_sort$ord
    }
    
    ord <- ord_list[[which.min(err_vec)]]
    
    if( !all(ord == seq(1,K))){
      new_t <- get_newtree(trace[[i]]$Ttree,ord)
      if(tree_check(new_t)){
        trace[[i]]$Ttree <- new_t
        trace[[i]]$theta <- trace[[i]]$theta[ord]
        trace[[i]]$phi <- trace[[i]]$phi[ord]
        trace[[i]]$L <- trace[[i]]$L[,ord]
        trace[[i]]$Z <- trace[[i]]$Z[,ord]
        trace[[i]]$Lo[,1] <- sapply(trace[[i]]$Lo[,1],function(x){ifelse(x==0,0,which(ord==x))})
        trace[[i]]$Zo[,1] <- sapply(trace[[i]]$Zo[,1],function(x){ifelse(x==0,0,which(ord==x))})
        trace[[i]]$C <- sapply(trace[[i]]$C,function(x){which(ord==x)})
      } else {
        Z_list[[length(Z_list)+1]] <- cur_Z 
      }
    }
  }
  
  return(trace)
}


# check if a vector represents a tree
tree_check <- function(tree)
{
  out <- T

  for(i in 1:length(tree)){
    if(i==1 && tree[1]!=0){return(F)}
    if(i==2 && tree[2]!=1){return(F)}
    if(tree[i]>=i){return(F)}
  }
  out
}

#############################################################
################ DISPLAY PARAMETERS ################
#############################################################
par_disp <- function(Params,MCMC_par)
{
  cat("Input information:\n")
  cat("Number of sequencing cells:",Params$N,"\n")
  cat("Number of loci:", Params$M,"\n")
  cat("Number of genome segments:", nrow(Params$segments),"\n")
  cat("\n")
  
  cat("Sampling parameters:\n")
  cat("Tuning samples:",MCMC_par$Ntune,"\n")
  cat("Burn-in samples:",MCMC_par$burnin,"\n")
  cat("Number of samples to be kept:", MCMC_par$Nsamp,"\n")
  cat("Number of MCMC chains:",MCMC_par$Nchain,
      ", with temperature increment:",MCMC_par$delta_T,"\n")
}



