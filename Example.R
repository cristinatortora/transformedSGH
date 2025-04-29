
### Data generated from a PTSGH
data1=rGHPT(300, 2,lambda = c(0.9,0.2))
data2=rGHPT(300, 2,lambda = c(0.5,0.2),mu=c(10,10))
data=(rbind(data1,data2))
plot(data)
l=c(rep(1,300),rep(2,300))

##Fitting MTGH and PTGH, it takes a few minutes
resM=MTGH(scale(data),2,initialization='km',iter.max=100,tol=0.001)
resP=PTGH(scale(data),2,initialization='km',iter.max=100,tol=0.001)

### Results
table(resM$labels,l)
table(resP$labels,l)


resM$BIC
resP$BIC

plot(data,col=resM$labels)
plot(data,col=resP$labels)

### Data generated from a MTSGH
data1=rGHMT(300, 2,lambda = c(0.9,0.2))
data2=rGHMT(300, 2,lambda = c(0.5,0.2),mu=c(10,10))
data=(rbind(data1,data2))
plot(data)
l=c(rep(1,300),rep(2,300))

##Fitting MTGH and PTGH, it takes a few minutes
resM=MTGH(scale(data),2,initialization='km',iter.max=100,tol=0.001)
resP=PTGH(scale(data),2,initialization='km',iter.max=100,tol=0.001)


table(resM$labels,l)
table(resP$labels,l)


resM$BIC
resP$BIC

plot(data,col=resM$labels)
plot(data,col=resP$labels)