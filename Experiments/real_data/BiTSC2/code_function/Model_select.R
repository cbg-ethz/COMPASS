####################################################################################
################ Bayesian Information Criterion ############################################################
################################################################################


BIC_fun <- function(X,D,onefolder){
  #### model selection
  rdata <- grep("K.*.Rdata",list.files(onefolder),value = T)
  #rdata <- list.files(onefolder)
  BIC <- list(K=c(),BIC=c())
  for(i in 1:length(rdata)){
    onerdata <- rdata[i]
    temp <- paste(onefolder,"/",onerdata,sep = "")
    load(temp)
    log_like <- mean(beta_log_like[1,])
    BIC$K <- c(BIC$K,Params$K)
    p <- 2*Params$K+5*Params$M+Params$N+5 #number of paramters: theta(K)+C(N)+pi(1)+Lo(M*2)+Zo(M*2)+g(M)+tree(K)+rho(1)+w(1)+s(1)+mu(1)
    n <- Params$N * Params$M *2  #number of samples
    bic <- -2*log_like+p*log(n)
    BIC$BIC <- c(BIC$BIC, bic)
  }
  plot(BIC$BIC ~ BIC$K, main = "BIC", xlab = "K", ylab = "BIC")
  save(BIC, file = paste(onefolder, "/BIC_model_selection.Rdata",sep = ""))
  BIC
}

