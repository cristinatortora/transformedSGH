library(MixGHD)
#source(MTGH)
#source(PTGH)
#source(rGH_Transformed)

#####GHD
data1=rGHD(300,2,alpha=c(2,-2))
data2=rGHD(300,2,alpha=c(-15,-10),mu=c(-3,-3))
data=(rbind(data1,data2))
plot(data)
l=c(rep(1,300),rep(2,300))
res=MTGH(scale(data),2,initialization='km',iter.max=100,tol=0.001)
resP=PTGH(scale(data),2,initialization='km',iter.max=100,tol=0.001)
resG=MGHD(data,G=2,max.iter=200,eps=0.001)

table(res$labels,l)
table(resP$labels,l)
table(resG@map,l)

res$BIC
resP$BIC
-resG@BIC



### PT SGH
data1=rGHPT(300, 2,lambda = c(0.9,0.2))
data2=rGHPT(300, 2,lambda = c(0.5,0.2),mu=c(10,10))
data=(rbind(data1,data2))
plot(data)
l=c(rep(1,300),rep(2,300))

res=MTGH(scale(data),2,initialization='km',iter.max=100,tol=0.001)
resP=PTGH(scale(data),2,initialization='km',iter.max=100,tol=0.001)
resG=MGHD(data,G=2,max.iter=200,eps=0.001)

table(res$labels,l)
table(resP$labels,l)
table(resG@map,l)

res$BIC
resP$BIC
-resG@BIC

### MT SGH
data1=rGHMT(300, 2,lambda = c(0.9,0.2))
data2=rGHMT(300, 2,lambda = c(0.5,0.2),mu=c(10,10))
data=(rbind(data1,data2))
plot(data)
l=c(rep(1,300),rep(2,300))

res=MTGH(scale(data),2,initialization='km',iter.max=100,tol=0.001)
resP=PTGH(scale(data),2,initialization='km',iter.max=100,tol=0.001)
resG=MGHD(data,G=2,max.iter=200,eps=0.001)

table(res$labels,l)
table(resP$labels,l)
table(resG@map,l)

res$BIC
resP$BIC
-resG@BIC