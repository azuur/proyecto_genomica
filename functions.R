require(here)
require(tidyverse)
require(magrittr)
require(infotheo)
require(minet)
require(GENIE3)
require(bnlearn)
require(tidyverse)
require(reshape2)
require(magrittr)
require(fastGeneMI)
require(genefilter)
require(lattice)
require(RColorBrewer)
require(latticeExtra)
require(vsn)
require(doParallel)
require(doRNG)
require(infotheo)
require(flare)
require(parallel)

#require(janitor)


## TRANSFORMATION FROM [0,infty) TO [0,1]
MIM_transform <- function(x){return(sqrt(1-exp(-2*x)))}


#TIGRESS ON ONE VARIABLE
#WE OUGHTA FIND THE SOURCE OF THIS. FOUND IT ONLINE
tigressIND<- function(y,x,alpha = 0.4,L = 2,R = 1000){
  n <- length(y);
  require(lars);
  #indexMat <- matrix(0,R,L)
  indexMat <- matrix(0,L,ncol(x))
  for(i in 1:floor(R/2)){
    indexVec <- sample(1:n,n);
    xr1 <- t(t(x[indexVec[1:floor(n/2)],])*runif(ncol(x),alpha,1));
    xr2 <- t(t(x[indexVec[(floor(n/2)+1):n],])*runif(ncol(x),alpha,1));
    result1 <- lars(x=xr1,y=y[indexVec[1:floor(n/2)]],type='lar',max.steps=L,use.Gram=FALSE)
    w1<-rev(order(result1$entry,decreasing=T)[1:L])
    indexMat[,w1][lower.tri(indexMat[,w1],diag=T)] <- indexMat[,w1][lower.tri(indexMat[,w1],diag=T)] + 1;
    result2 <- lars(x=xr2,y=y[indexVec[(floor(n/2)+1):n]],type='lar',max.steps=L,use.Gram=FALSE)
    w2<-rev(order(result2$entry,decreasing=T)[1:L])
    indexMat[,w2][lower.tri(indexMat[,w2],diag=T)] <- indexMat[,w2][lower.tri(indexMat[,w2],diag=T)] + 1;
  }
  return(colMeans(indexMat/R))
}

##TIGRESS ON DATASET
tigress<-function(x,alpha = 0.4,L = 2,R = 1000, verbose=FALSE){ #optimal perfomance has R = 8000
  
  #OJOOO
  x <- scale(x)
  
  
  nvar <- ncol(x)
  net <- matrix(0,ncol=nvar,nrow=nvar)
  for(i in 1:nvar){
    if(verbose==T){print(sprintf("Bootstrapping least-angle regression on variable %s of %s.",i,nvar))}
    net[-i,i] <- tigressIND(c(x[,i]),x[,-i],alpha = alpha,L = L,R = R)
  }
  return(net)
}

##NARROMI
RO_alg <- function(i,x,MIM,lambda,eps){
  candidate_TF <- unname(which(MIM[,i]!=0))
  zeroes <- 1
  #    huh <-0 
  while(zeroes > 0){
    if(length(candidate_TF)==0){break}
    RO <- slim(X = as.matrix(x[,candidate_TF]), Y = x[,i], 
               method = "lq", q=1, lambda=lambda,verbose = F, nlambda = 1)#,...)
    zero <- abs(RO$beta) < eps
    candidate_TF <- candidate_TF[!zero]
    zeroes <- sum(zero)
    #      huh <- huh+1
    #      print(huh)
  }
  res <- list(candidate_TF=candidate_TF,beta=RO$beta)
  return(res)
}

narromi<-function(x, MIM =NULL, t = 0.6, estimator = "mi.empirical", disc = "equalfreq",lambda = 1,theta = 0.05, eps = 0.05,verbose = FALSE, parallel = F, numcores = NULL,...){

  if(is.null(MIM)){MIM <- build.mim(x,estimator,...)}

  #VARIABLE DEF
  ##OJOOO!! ESTÁ ESCALADO
  x <- scale(x)
  x <- as.matrix(x)
  nvar <- nrow(MIM)
  if(theta == "choose"){theta <- quantile(sort(abs(MIM))[-c(1:nvar)],0.15)}
  if(verbose == TRUE){ print("Building/finding MI matrix.")}
  MIM <- as.matrix(MIM)
  MIM[abs(MIM)<theta] <- 0
  diag(MIM) <- 0
  net <- matrix(0,nrow = nvar, ncol = nvar)
  
  RO_on_var_num <- function(i){
    if(verbose == TRUE){print( sprintf("Executing RO step on variable %s of %s.",i,nvar) )}
    vect <- rep(0,nvar)
    RO <- RO_alg(i,x=x,MIM=MIM,lambda=lambda,eps=eps)
    candidate_TF <- RO$candidate_TF
    beta <- RO$beta
    if(length(candidate_TF)>0){
      transformed_MI <- MIM_transform(MIM[,i][candidate_TF])
      vect[candidate_TF]<-sign(RO$beta)*(t*transformed_MI + (1-t)*abs(RO$beta))
    }
    return(vect)
  }
#  684486094
  if(parallel == T){
    if(is.null(numcores)){ numcores <- parallel::detectCores() - 1 }
    cl <- makeCluster(numcores)
    parallel::clusterExport(cl = cl,list("RO_alg","slim","MIM_transform"))
    net <- parSapply(cl = cl,X = 1:nvar, FUN = RO_on_var_num)
    stopCluster(cl)
  } else{
    net <- sapply(X = 1:nvar, FUN = RO_on_var_num)
  }

  dimnames(net) <- dimnames(MIM)
  return(net)
}



##DEPRECATED IN FAVOR OF PARALLELIZED VERSION
# narromi1<-function(x, MIM =NULL, t = 0.6, estimator = "mi.empirical", disc = "equalfreq",lambda = 1,theta = 0.05, eps = 0.05,verbose = FALSE, ...){
#   require("minet")
#   require("infotheo")
#   require("flare")
#   
#   if(is.null(MIM)){MIM <- build.mim(x,estimator,...)}
#   
#   ##OJOOO!!
#   x <- scale(x)
#   
#   nvar <- nrow(MIM)
#   
#   if(theta == "choose"){theta <- quantile(sort(abs(MIM))[-c(1:nvar)],0.15)}
#   
#   if(verbose == TRUE){ print("Building/finding MI matrix.")}
#   MIM[abs(MIM)<theta] <- 0
#   diag(MIM) <- 0
#   net <- matrix(0,nrow = nvar, ncol = nvar)
#   
#   for(i in 1:nvar){
#     if(verbose == TRUE){print( sprintf("Executing RO step on variable %s of %s.",i,nvar) )}
#     
#     candidate_TF <- unname(which(MIM[,i]!=0))
#     zeroes <- 1
#     #    huh <-0
#     while(zeroes > 0){
#       if(length(candidate_TF)==0){break}
#       RO <- slim(X = as.matrix(x[,candidate_TF]), Y = x[,i],
#                  method = "lq", q=1, lambda=lambda,verbose = F, nlambda = 1)#,...)
#       zero <- abs(RO$beta) < eps
#       candidate_TF <- candidate_TF[!zero]
#       zeroes <- sum(zero)
#       #      huh <- huh+1
#       #      print(huh)
#     }
#     
#     if(length(candidate_TF)>0){
#       transformed_MI <- MIM_transform(MIM[,i][candidate_TF])
#       net[,i][candidate_TF]<-sign(RO$beta)*(t*transformed_MI + (1-t)*abs(RO$beta))
#     }
#     
#     
#   }
#   return(net)
# }
# 



##SAVE AND LOAD R OBJECTS
saveobj<-function(x, gitignore = F){
  dir <- "Robjects/"
  if(gitignore == T){dir <- "Robjects/Results"}
  assign(x,get(x,envir = .GlobalEnv),envir = environment())
  saveRDS(get(x,envir = environment()),here(paste0(dir,x,".RDS")))
} 
readobj<-function(x) assign(x,readRDS(here(paste0("/Robjects/",x,".RDS"))),envir = .GlobalEnv)

savelist <- function(...){
  list <- c(...)
  for(i in list) saveobj(i)
}
readlist <- function(...){
  list <- c(...)
  for(i in list) readobj(i)
}

robject_list <- function(){
  l <- list.files(here("robjects"))
  l_nonames <- tools::file_path_sans_ext(l)
  return(l_nonames)
} 

