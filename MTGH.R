


################################
### MAIN
##############################

#########################################################################
#                         Model based clustering using                  #
#  the Manly transformed symmetric generalized hyperbolic distribution  #
#########################################################################




MTGH<-function(data,G=2,iter.max=100,tol=0.001,initialization="pam"){
  ##data              numerical dataset
  ##G                 number of clusters
  ##iter.max          maximum number of iterations
  ##tol               tolerance for Aitken convergence
  ##initialization    initialization method, options: "mclust" Gaussian mixture models, "km" k-means, 
  ##                    "pam"k-medoids (default)
  
  X=as.matrix(data)
  n=nrow(X)
  p=ncol(X)
  
  ####initialization
  if(G>1){
    # if(initialization=="random.soft"){
    #   
    #   z  <- array(runif(n*k),c(n,k)) # soft posterior probabilities (no-normalized) (n x k)
    #   z  <- z/rowSums(z)             # soft posterior probabilities (n x k)
    #   
    # }
    # 
    # if(initialization=="random.hard"){
    #   
    #   z  <- t(rmultinom(n, size = 1, prob=rep(1/k,k)))  # hard posterior probabilities (n x k)
    #   
    # }
    # 
    # if(initialization=="manual"){ # z.start can be both soft and hard initialization
    #   
    #   z  <- start.z # soft or hard classification
    #   
    # }
    
    if(initialization=="mclust"){
      
      if(d==1)
        mclustfit <- mclust::Mclust(data=X, G=G, modelNames="V")
      if(d>1)
        mclustfit <- mclust::Mclust(data=X, G=G, modelNames="VVV")
      
      z <- mclustfit$z
      mu=mclustfit$parameters$mean
      
    }
    
    if(initialization=="km"){
      
      km       <- kmeans(X,G)
      groups   <- km$cluster
      z        <- mclust::unmap(classification=groups)
      
      
      mu=km$centers
    }
    
    if(initialization=="pam"){
      
      lk     <- cluster::pam(x=X,k=G)
      groups <- lk$clustering
      z      <- mclust::unmap(classification=groups)
      mu=lk$medoids
    }
    
    
    
    
  }
  if(G==1){
    groups   <- rep(1,n)
    z        <- matrix(groups,nrow=n,ncol=k)
    mu       <-apply(X,2,mean)
  }
  
  prior    <-apply(z,2,mean)
  ss    <-diag(p)
  sigma <-array(ss,c(p,p,G))
  cpl   <- matrix(c(1,-1/2),G,2,1)
  PTl   <- matrix(0.9,G,p)
  
  y <- array(0,c(n,p,G))
  
  llik=NULL
  for(iter in 1:3){
    
    ##### z update
    
    denZ=dMTSGH(X,G,PTl=PTl,mu=mu,sigma=sigma,omega=cpl[,1],lambda=cpl[,2],prior=prior)
    z=matrix(0,n,G)
    for (k in 1:G){
      ss=array(sigma[,,k],c(p,p,1))
      numZ=prior[k]*dMTSGH(X,G=1,PTl=matrix(PTl[k,],nrow=1,ncol=p),mu=matrix(mu[k,],1,p),sigma=ss,omega=matrix(cpl[k,1],1,1),lambda=matrix(cpl[k,2],1,1),prior=1)
      z[,k]=numZ/denZ
    }
    
    
    
    
    for(k in 1:G){
      #Data Transformation
      
      y[,,k]=Manly(X=X,PTl[k,])
      
      #update mu sigma cpl
      par=list()
      par$mu=mu[k,]
      par$sigma=sigma[,,k]
      par$cpl=cpl[k,]
      par = updateEM(x= y[,,k], par=par, weights=z[,k], invS=NULL)
      mu[k,]=par$mu
      sigma[,,k]=par$sigma
      cpl[k,]=par$cpl
    }
    
    #update pi
    prior = apply(z,2,mean)
    
    #update PTl
    for(k in 1:G){
      res = optim(par=PTl[k,], fn=fobjPTl_MT, data=X,mu=as.vector(mu[k,]),sigma=sigma[,,k],cpl=cpl[k,],z=z[,k],method = "Nelder-Mead",
                  control = list(fnscale=-1))
      PTl[k,] = res$par
    }
    
    llik=c(llik,sum(log(dMTSGH(X,G,PTl=PTl,mu=mu,sigma=sigma,omega=cpl[,1],lambda=cpl[,2],prior=prior))))
  }
  
  
  
  iteration=3
  check=0
  
  while (check==0) {
    iteration=iteration+1
    
    #  if(iteration %% 10 == 0) cat("Iteration",iteration,"\n")
    
    ##### z update
    
    denZ=dMTSGH(X,G,PTl=PTl,mu=mu,sigma=sigma,omega=cpl[,1],lambda=cpl[,2],prior=prior)
    z=matrix(0,n,G)
    for (k in 1:G){
      ss=array(sigma[,,k],c(p,p,1))
      numZ=prior[k]*dMTSGH(X,G=1,PTl=matrix(PTl[k,],nrow=1,ncol=p),mu=matrix(mu[k,],1,p),sigma=ss,omega=matrix(cpl[k,1],1,1),lambda=matrix(cpl[k,2],1,1),prior=1)
      z[,k]=numZ/denZ
      # }
      
      
      #for(k in 1:G){
      #Data Transformation
      
      y[,,k]=Manly(X=X,PTl[k,])
      
      par=list()
      par$mu=mu[k,]
      par$sigma=sigma[,,k]
      par$cpl=cpl[k,]
      par = updateEM(x= y[,,k], par=par, weights=z[,k], invS=NULL)
      mu[k,]=par$mu
      sigma[,,k]=par$sigma
      cpl[k,]=par$cpl
      
      res = optim(par=PTl[k,], fn=fobjPTl_MT, data=X,mu=as.vector(mu[k,]),sigma=sigma[,,k],cpl=cpl[k,],z=z[,k],method = "Nelder-Mead",
                  control = list(fnscale=-1))
      PTl[k,] = res$par
    }
    prior = apply(z,2,mean)
    llik=c(llik,sum(log(dMTSGH(X,G,PTl=PTl,mu=mu,sigma=sigma,omega=cpl[,1],lambda=cpl[,2],prior=prior))))
    
    if(iteration == iter.max |  getall(llik)==TRUE ) check <- 1#
  }
  # print(all(llik==sort(llik)))
  labels=apply(z,1,which.max)
  # pairs(data,col=labels)
  
  
  
  
  # -------------------- #
  # Number of parameters #
  # -------------------- #
  final.loglik = llik[iteration ]
  
  npar <- (G-1) + p*G + p*G + p*(p+1)/2*G + 2*G 
  
  # -------------------- #
  # Information criteria #
  # -------------------- #
  
  AIC <- - 2*final.loglik + 2*npar
  BIC <- - 2*final.loglik + npar*log(n)
  
  KIC <- - 2*final.loglik + 3*(npar+1)
  # KICc <- - 2*final.loglik + 2*(npar + 1)*n/(n-npar -2) - n*digamma((n-npar)/2) + n*log(n/2)
  
  AIC3 <- - 2*final.loglik + 3*npar
  CAIC <- - 2*final.loglik + npar*(1+log(n))
  AICc <- - 2*final.loglik + 2*npar*n/(n - npar -1)
  
  ent <- apply(z,1,max)
  ICL <- BIC - sum(ent*log(ent))
  
  AWE <- - 2*(final.loglik + sum(ent*log(ent))) + 2*npar*(3/2 + log(n))
  CLC <- - 2*final.loglik + 2*sum(ent*log(ent))
  
  
  return(list(llik=llik,mu=mu,sigm=sigma, PTl=PTl,cpl=cpl,prior=prior,z=z,
              labels=labels,
              AIC            = AIC,
              BIC            = BIC,
              ICL            = ICL,
              KIC            = KIC,
              #   KICc           = KICc,
              AWE            = AWE,
              AIC3           = AIC3,
              CAIC           = CAIC,
              AICc           = AICc,
              #AICu           = AICu,
              CLC            = CLC
  ))
}


##########################
## Manly Transformation ##
##########################

Manly <- function(X,lambda = 0.1){
  
  n <- length(X)
  Y <- X
  lX=sweep(X,2,FUN='*',lambda)
  #if(lambda != 0){
  Yt=(exp(lX)-1)
  Y=sweep(Yt,2,FUN='/',lambda)
  #  Y=(exp(lX)-1)/lambda
  #}library(cranlogs)

  
  return(Y)
  
}

Manly2 <- function(X,lambda = 0.1){
  
  n <- length(X)
  Y <- X
  J <- rep(1,n)
  lX=sweep(X,2,FUN='*',lambda)
  #if(lambda != 0){
  Yt=(exp(lX)-1)
  Y=sweep(Yt,2,FUN='/',lambda)
  J=exp(lX)
  # if(lambda != 0){
  #   Y=(exp(lambda*X)-1)/lambda
  #   J=exp(lambda*X)
  # }
  
  return(list(Y=Y,J=J))
  
}


###############################################
#####. density of a mixture of Manly transformed SGH
###############################################
dMTSGH<-function(data,G,PTl,mu,sigma,omega,lambda,prior){
  p=ncol(data)
  n=nrow(data)
  
  dens  <- array(1,c(n,G),dimnames=list(1:n,paste("comp.",1:G,sep=""))) 
  ymat=matrix(0,n,p)
  jacmat=matrix(0,n,p)
  for(k in 1:G){
    jac=Manly2(data,PTl[k,])
    ymat=jac$Y
    jacmat=jac$J
    # for(j in 1:p){
    #   jac=Manly2(data[,j],PTl[k,j])
    #   ymat[,j]=jac$Y
    #   jacmat[,j]=jac$J
    # }
    dens[,k]<-dSGH(matrix(ymat,n,p),p,matrix(mu[k,],1,p),sigma[,,k], omega[k],lambda[k])*apply(jacmat,1,prod)*  prior[k]
  }
  if(k>1){dens=apply(dens,1,sum)}
  PDF=dens
  PDF <- (PDF<10^(-323))*10^(-323)+(PDF>10^(-323))*PDF
  return(PDF)
}

#########################
### PTl update
########################



fobjPTl_MT = function(par,data,mu,cpl,sigma,z){
  ll=dMTGH1k(data=data,PTl=par,mu=mu,sigma=sigma,cpl=cpl,z=z)
  ll=ll[which(ll!=0)]
  ll = sum(log(ll))
  #ll = prod((dPmscn1k(data=data,k=k,PTl=par,mu=mu,alpha=alpha,eta=eta,Lambda=Lambda,Gamma=Gamma,z=z)))
  return(ll)
}




dMTGH1k<-function(data,PTl,mu,sigma,cpl,z){
  p=ncol(data)
  n=nrow(data)
  
  dens  <- rep(1,n)
  ymat=matrix(0,n,p)
  jacmat=matrix(0,n,p)
  omega=cpl[1]
  lambda=cpl[2]

    jac=Manly2(data,PTl)
    ymat=jac$Y
    jacmat=jac$J

  dens<-dSGH(ymat,p,mu,sigma, omega,lambda)*apply(jacmat,1,prod)
  
  # dens <- (dens<10^(-323))*10^(-323)+(dens>10^(-323))*dens
  dens=(dens)^z#exp(log(dens)*z)
  
  #dens <- (dens<10^(-323))*10^(-323)+(dens>10^(-323))*dens
  return(dens)
}


