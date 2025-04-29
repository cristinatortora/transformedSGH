#########################################################################
#                 pseudo random number generation using                 #
# rGHPT: power transformed symmetric generlized hyperbolic distribution #
# rGHMT: Manly transformed symmetric generlized hyperbolic distribution #
#########################################################################

library(MixGHD)

rGHMT<- function(n, d,lambda = rep(0.5,d), mu = rep(10,d), Sigma=diag(d),omega=1,LA=0.5){
  
  gen=rGHD(n,p=d, mu=mu,sigma=Sigma,omega=omega,lambda=LA)
  
  Y=gen 
  p=ncol(Y)
  X=Y
  for(j in 1:p){
    if(lambda[j]!=0)
      X[,j]=lambda[j]^(-1)*log(Y[,j]*lambda[j]+1)
  }
  return(X=X)
  
}


rGHPT<- function(n, d,lambda = rep(0.5,d), mu = rep(10,d), Sigma=diag(d),omega=1,LA=0.5){
  

  gen=rGHD(n,p=d, mu=mu,sigma=Sigma,omega=omega,lambda=LA)
  
  Y=gen 
  p=ncol(Y)
  X=matrix(0,n,p)
  for(j in 1:p){
    X[,j]=PTinv(Y[,j],lambda[j])
  }
  return(X=X)
  
}

PTinv <- function(Y,lambda = 1){
  
  n <- length(Y)
  X <- numeric(n)
  
  for(i in 1:n){
    if(Y[i] >= 0 & lambda != 0)
      X[i] <-(Y[i]*lambda+2)^(1/lambda)-1
    if(Y[i] >= 0 & lambda == 0)
      X[i] <- exp(Y[i]) - 1
    if(Y[i] < 0 & lambda != 2)
      X[i] <- 1-(1-Y[i]*(2-lambda))^(1/(2-lambda))
    if(Y[i] < 0 & lambda == 2)
      X[i] <- 1- exp(-Y[i])
  }
  
  return(X)
  
}