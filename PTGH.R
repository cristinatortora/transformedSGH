#########################################################################
#                         Model based clustering using                  #
#  the power transformed symmetric generalized hyperbolic distribution  #
#########################################################################
library(MASS)
library(Bessel)
library(ghyp)


################################
### MAIN
##############################

PTGH<-function(data,G=2,iter.max=100,tol=0.0001,initialization="pam"){
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
  PTl   <- matrix(1,G,p)
  
  y <- array(0,c(n,p,G))
  
  llik=NULL
  for(iter in 1:3){
    
    ##### z update
    
    denZ=dTSGH(X,G,PTl=PTl,mu=mu,sigma=sigma,omega=cpl[,1],lambda=cpl[,2],prior=prior)
    z=matrix(0,n,G)
    for (k in 1:G){
      ss=array(sigma[,,k],c(p,p,1))
      numZ=prior[k]*dTSGH(X,G=1,PTl=matrix(PTl[k,],nrow=1,ncol=p),mu=matrix(mu[k,],1,p),sigma=ss,omega=matrix(cpl[k,1],1,1),lambda=matrix(cpl[k,2],1,1),prior=1)
      z[,k]=numZ/denZ
    }
    
    
    
    
    for(k in 1:G){
      #Data Transformation
      for( h in 1:p){
        y[,h,k]=PT(X=X[,h],PTl[k,h])
      }
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
      res = optim(par=PTl[k,], fn=fobjPTl, data=X,mu=as.vector(mu[k,]),sigma=sigma[,,k],cpl=cpl[k,],z=z[,k],method = "Nelder-Mead",
                  control = list(fnscale=-1))
      PTl[k,] = res$par
    }
    
    llik=c(llik,sum(log(dTSGH(X,G,PTl=PTl,mu=mu,sigma=sigma,omega=cpl[,1],lambda=cpl[,2],prior=prior))))
  }
  
  
  
  iteration=3
  check=0
  
  while (check==0) {
    iteration=iteration+1
    
    #if(iteration %% 10 == 0) cat("Iteration",iteration,"\n")
    
    ##### z update
    
    denZ=dTSGH(X,G,PTl=PTl,mu=mu,sigma=sigma,omega=cpl[,1],lambda=cpl[,2],prior=prior)
    z=matrix(0,n,G)
    for (k in 1:G){
      ss=array(sigma[,,k],c(p,p,1))
      numZ=prior[k]*dTSGH(X,G=1,PTl=matrix(PTl[k,],nrow=1,ncol=p),mu=matrix(mu[k,],1,p),sigma=ss,omega=matrix(cpl[k,1],1,1),lambda=matrix(cpl[k,2],1,1),prior=1)
      z[,k]=numZ/denZ
    }
    
    
    for(k in 1:G){
      #Data Transformation
      for( h in 1:p){
        y[,h,k]=PT(X=X[,h],PTl[k,h])
      }
      par=list()
      par$mu=mu[k,]
      par$sigma=sigma[,,k]
      par$cpl=cpl[k,]
      par = updateEM(x= y[,,k], par=par, weights=z[,k], invS=NULL)
      mu[k,]=par$mu
      sigma[,,k]=par$sigma
      cpl[k,]=par$cpl
    }
    prior = apply(z,2,mean)
    
    
    for(k in 1:G){
      res = optim(par=PTl[k,], fn=fobjPTl, data=X,mu=as.vector(mu[k,]),sigma=sigma[,,k],cpl=cpl[k,],z=z[,k],method = "Nelder-Mead",
                  control = list(fnscale=-1))
      PTl[k,] = res$par
    }
    
    llik=c(llik,sum(log(dTSGH(X,G,PTl=PTl,mu=mu,sigma=sigma,omega=cpl[,1],lambda=cpl[,2],prior=prior))))
    
    if(iteration == iter.max |  getall(llik)==TRUE) check <- 1
  }
  #print(all(llik==sort(llik)))
  labels=apply(z,1,which.max)
  #pairs(data,col=labels)
  
  
  
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
  KICc <- - 2*final.loglik + 2*(npar + 1)*n/(n-npar -2) - n*digamma((n-npar)/2) + n*log(n/2)
  
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
              KICc           = KICc,
              AWE            = AWE,
              AIC3           = AIC3,
              CAIC           = CAIC,
              AICc           = AICc,
              #AICu           = AICu,
              CLC            = CLC
  ))
}

##########################
## Power Transformation ##
##########################

PT <- function(X,lambda = 1){
  
  n <- length(X)
  Y <- numeric(n)
  
  for(i in 1:n){
    if(X[i] >= 0 & lambda != 0){
      Y[i] <- 1/lambda*((X[i] + 1)^(lambda) - 1)}
#    if(X[i] >= 0 & lambda == 0)
#      Y[i] <- log(X[i] + 1)
   # if(X[i] < 0 & lambda != 2)
    else{
      Y[i] <- 1/(2 - lambda)*(1 - (1 - X[i])^(2 - lambda))}
 #   if(X[i] < 0 & lambda == 2)
  #    Y[i] <- -log(1-X[i])
  }
  
  return(Y)
  
}

PT2 <- function(X,lambda = 1){
  
  n <- length(X)
  Y <- numeric(n)
  J <- numeric(n)
  
  for(i in 1:n){
    if(X[i] >= 0 & lambda != 0){
      Y[i] <- 1/lambda*((X[i] + 1)^(lambda) - 1)
      J[i] <- abs((X[i]+1)^(lambda-1))
    }
    # if(X[i] >= 0 & lambda == 0){
    #   Y[i] <- log(X[i] + 1)
    #   J[i] <- abs((X[i]+1)^(-1))
    # }
   # if(X[i] < 0 & lambda != 2){
      else{
      Y[i] <- 1/(2 - lambda)*(1 - (1 - X[i])^(2 - lambda))
      J[i] <- abs((2-lambda)^(-1)*((2-lambda)*(1-X[i])^(1-lambda)))
    }
    # if(X[i] < 0 & lambda == 2){
    #   Y[i] <- -log(1-X[i])
    #   J[i] <- abs((1-X[i])^(-1))
    # }
  }
  
  return(list(Y=Y,J=J))
  
}

###############################
##### density of a single SGH
###############################
dSGH <- function(data, p, mu=rep(0,p),sigma=diag(p),omega=1,lambda=0.5, log=FALSE) {
  # x  is a n * p matrix
  # generalized hyperbolic distribution
  if(!is.matrix(data))data=matrix(data,length(data)/p,p)
  x=data
  
  # mu    = par$mu
  # sigma = par$sigma
  # alpha = par$alpha
  
  #check for input error
   d = length(mu)
  omega = exp(log(omega))
  
  
  
  #distances between points and mu
  invS = ginv(sigma)
  pa = omega
  mx = omega + mahalanobis(x, center=mu, cov=invS, inverted=TRUE)
  pa2=rep(pa,length(mx))
  kx = sqrt(mx*pa2)
  xmu = sweep(x,2,mu)
  
  lvx = matrix(0, nrow=nrow(x), 2)
  lvx[,1] = (lambda - d/2)*log(kx)
  lvx[,2] = log(besselK( kx, nu=lambda-d/2, expon.scaled =TRUE)) - kx

  
  
  lv = numeric(6)
  if(is.nan(log(det( sigma )))){sigma =diag(ncol(mu))}
  
  lv[1] = -1/2*log(det( sigma )) -d/2*(log(2)+log(pi))
  lv[2] =  omega - log(besselK( omega, nu=lambda, expon.scaled =TRUE))
  lv[3] = -lambda/2*( log(1) )
  lv[4] = lambda*log(omega) *0
  lv[5] = (d/2 - lambda)*log( pa )
  
  val = apply(lvx,1,sum) + sum(lv)
  if (!log) val = exp( val )
  
  return(val)
}


###############################################
#####. density of a mixture of transformed SGH
###############################################
dTSGH<-function(data,G,PTl,mu,sigma,omega,lambda,prior){
  p=ncol(data)
  n=nrow(data)
  
  dens  <- array(1,c(n,G),dimnames=list(1:n,paste("comp.",1:G,sep=""))) 
  ymat=matrix(0,n,p)
  jacmat=matrix(0,n,p)
  for(k in 1:G){
    for(j in 1:p){
      jac=PT2(data[,j],PTl[k,j])
      ymat[,j]=jac$Y
      jacmat[,j]=jac$J
    }
   dens[,k]<-dSGH(ymat,p,mu[k,],sigma[,,k], omega[k],lambda[k])*apply(jacmat,1,prod)*  prior[k]
  }
  if(k>1){dens=apply(dens,1,sum)}
  PDF=dens
  PDF <- (PDF<10^(-323))*10^(-323)+(PDF>10^(-323))*PDF
  return(PDF)
}


###############################################
##### EM SGH
###############################################

updateEM <- function(x, par, weights=NULL, invS=NULL) {
  ####computation of mu and alpha sigma cpl
  ##########piu importante intervenire qui!!
  if (is.null(weights)) weights=rep(1,nrow(x))
  if (is.null(invS)) invS=ginv(par$sigma)
  
  # expectations of w, 1/w, log(w) given x
  abc = gig2SGH(x=x, par=par, invS=invS)
  d = length(par$mu)
  
  sumw = sum(weights)
  ABC = apply(abc,2,weighted.sum, wt=weights)/sumw

    A = ABC[1]
    B = ABC[2]
    u = (B - abc[,2])*weights
    t = (A*abc[,2]-1)*weights
   
    
    mu.new   = apply(x, 2, weighted.mean, w=abc[,2]*weights) 


  
  R = cov.wt(x, wt=abc[,2]*weights, center=mu.new, method="ML")$cov*ABC[2] #returns a list containing weighted covariance matrix and the mean of the data

  for(i in 1:ncol(R)){if(R[i,i]<0.00001){R[i,i]=0.00001}}
  par.old = c(  log(par$cpl[1]) , par$cpl[2])

  par.ol = c(exp(par.old[1]), par.old[2])
  a = updateol(ol=par.ol, ABC=ABC, n=2)

  cpl.new =  c( a[1], a[2])
  new.par = list(mu=mu.new, sigma=R, cpl=cpl.new )
  return(new.par)
}


###############################################
###### E step SGH
###############################################


gig2SGH <- function(x=NULL, par=NULL, invS=NULL) {
  #  returns a matrix with dim length(a) x 3
  d = length(par$mu)
  
  #if (is.null(invS)) invS = ginv(par$sigma)
  omega = exp(  log(par$cpl[1])  )
  
  a1 = omega 
  b1 = omega + as.numeric(mahalanobis(x, center=par$mu, cov=invS, inverted=TRUE))
  v1 = par$cpl[2]-d/2
  
  val = gigGH(b=b1,a=a1, v=v1)
  return(val)
}

gigGH <- function(a=NULL,b=NULL,v=NULL) {
  # returns a matrix with dim length(a) x 3 stima le yi
  sab =  sqrt(a*b)
  kv1 = besselK( sab, nu=v+1, expon.scaled =TRUE)
  kv  = besselK( sab, nu=v, expon.scaled =TRUE)
  kv12 = kv1/kv
  
  sb.a = sqrt(b/a)
  w    = kv12*sb.a
  invw = kv12*1/sb.a - 2*v/b
  logw = log(sb.a) + numDeriv::grad( logbesselKv, x=rep(v,length(sab)), y=sab, method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
  
  
  val =cbind(w,invw,logw)
  return(val)
}

weighted.sum <- function(z,wt,...) return( sum(z*wt,...) )

###Bessel function
logbesselKv <- function(x, y) { log(besselK(x=y, nu=x, expon.scaled=TRUE)) - log(y)}



###########################
### cpl update
############################
updateol <- function(ol=NULL, ABC=NULL, n=1) {
  
  for (i in 1:n) {
    if (ABC[3] == 0) {
      ol[2] = 0
    } else {
      # if(ol[1]<=0){ol[1]=eps}
      bv = grad( logbesselKv, x=ol[2], y=ol[1], method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      ol[2] = ABC[3]*(ol[2]/bv)
    }
    
    
    lam0 = ol[2]
    omg0 = ol[1]
    Rp = RlamU(omg0,lam=+lam0)
    Rn = RlamU(omg0,lam=-lam0)
    f1 = Rp + Rn - (ABC[1]+ABC[2])
    f2 = ( Rp^2 - (2*lam0+1)/omg0 *Rp -1 ) + ( Rn^2 - (2*(-1*lam0)+1)/omg0 *Rn -1 )
    # note, it is suppose to be f1/2 and f2/2
    if ( ol[1] - f1/f2 > 0 ) ol[1] = ol[1] - f1/f2
  }
  
  return(ol)
}
RlamU <- function(x, lam=NULL) {
  
  v1 = besselK(x, nu= lam+1, expon.scaled=TRUE)
  v0 = besselK(x, nu= lam, expon.scaled=TRUE)
  val = v1/v0
  return(val)
}

#########################
### PTl update
########################



fobjPTl = function(par,data,mu,cpl,sigma,z){
  ll=dTGH1k(data=data,PTl=par,mu=mu,sigma=sigma,cpl=cpl,z=z)
  ll=ll[which(ll!=0)]
  ll = sum(log(ll))
  #ll = prod((dPmscn1k(data=data,k=k,PTl=par,mu=mu,alpha=alpha,eta=eta,Lambda=Lambda,Gamma=Gamma,z=z)))
  return(ll)
}




dTGH1k<-function(data,PTl,mu,sigma,cpl,z){
  p=ncol(data)
  n=nrow(data)
  
  dens  <- rep(1,n)
  ymat=matrix(0,n,p)
  jacmat=matrix(0,n,p)
  omega=cpl[1]
  lambda=cpl[2]
  for(h in 1:p){
    jac=PT2(data[,h],PTl[h])
    ymat[,h]=jac$Y
    jacmat[,h]=jac$J
  }
  dens<-dSGH(ymat,p,mu,sigma, omega,lambda)*apply(jacmat,1,prod)

  # dens <- (dens<10^(-323))*10^(-323)+(dens>10^(-323))*dens
  dens=(dens)^z#exp(log(dens)*z)
  
  #dens <- (dens<10^(-323))*10^(-323)+(dens>10^(-323))*dens
  return(dens)
}

#####Stopping criteria

getall <- function(loglik) {
  if (length(loglik) <3) stop("must have at least 3 likelihood values")
  n = length(loglik)
  lm1 = loglik[n]
  lm  = loglik[(n-1)]
  lm_1  = loglik[(n-2)]
  am = (lm1 - lm)/(lm - lm_1)
  lm1.Inf = lm + (lm1 - lm)/(1-am)
  val = lm1.Inf - lm
  if (is.nan(val)) val=0
  if (val < 0) val= 1
  return( val )
}

