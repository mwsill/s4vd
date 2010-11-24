# TODO: Add comment
# 
# Author: martin
###############################################################################
rm(list=ls())
source("~/workspace/s4vd/pkg/R/aaa-s4vd.R")
source("~/workspace/s4vd/pkg/R/s4vd2.0.R")
source("~/workspace/s4vd/pkg/R/stabpathplot.R")
source("~/workspace/s4vdsim/datagen.R")
source("~/workspace/s4vd/pkg/R/ssvdBC.R")
load("~/workspace/s4vdsim/epenHD.RData")
source("~/workspace/s4vdsim/heatmap.r")


#ependymoma data
seed <- 24122010
X <- t(scale(t(HDsel)))
set.seed(seed)
res5 <- s4vd(X,pcer=0.5,steps=500,c(0.6,0.65),nbiclust=10,merr=0.05,pointwise=TRUE,gamm=0,col.overlap=FALSE)
set.seed(seed)
res4 <- s4vd(X,pcer=0.4,steps=500,c(0.6,0.65),nbiclust=10,merr=0.05,pointwise=TRUE,gamm=0,col.overlap=FALSE)
set.seed(seed)
res3 <- s4vd(X,pcer=0.3,steps=500,c(0.6,0.65),nbiclust=10,merr=0.05,pointwise=TRUE,gamm=0,col.overlap=FALSE)
set.seed(seed)
res2 <- s4vd(X,pcer=0.2,steps=500,c(0.6,0.65),nbiclust=10,merr=0.05,pointwise=TRUE,gamm=0,col.overlap=FALSE)
set.seed(seed)
res25 <- s4vd(X,pcer=0.25,steps=500,c(0.6,0.65),nbiclust=10,merr=0.05,pointwise=TRUE,gamm=0,col.overlap=FALSE)
set.seed(seed)
res1 <- s4vd(X,pcer=0.1,steps=500,c(0.6,0.65),nbiclust=10,merr=0.05,pointwise=TRUE,gamm=0,col.overlap=FALSE)
set.seed(seed)
res15 <- s4vd(X,pcer=0.15,steps=500,c(0.6,0.65),nbiclust=10,merr=0.05,pointwise=TRUE,gamm=0,col.overlap=FALSE)



library(marray)
cols <-  maPalette(low = "blue", high = "red", mid="white", k =50)
heatmapBC(X,res1,col=cols,order=TRUE)
heatmapBC(X,res15,col=cols,order=TRUE)
heatmapBC(X,res2,col=cols,order=TRUE)
heatmapBC(X,res25,col=cols,order=TRUE)
heatmapBC(X,res3,col=cols,order=TRUE)
heatmapBC(X,res4,col=cols,order=TRUE)
heatmapBC(X,res5,col=cols,order=TRUE)

set.seed(12345)
resnp <- s4vd(X,pcer=0.1,steps=300,nbiclust=10,merr=0.05,pointwise=FALSE,gamm=0)



thr2 <- ((vc[[4]]^2/((n*0.25)*n))+1)/2
orange <- which(apply(vc[[2]][,1:sum(thr2<thr)],1,max)>thr)
thr3 <- ((vc[[4]]^2/((n*0.5)*n))+1)/2
green <- which(apply(vc[[2]][,1:sum(thr3<thr)],1,max)>thr)
cols <- rep("black",nrow(vc[[2]]))
cols[red] <- "red"
cols[orange] <- ifelse(cols[orange]=="black","orange",cols[orange])
cols[green] <- ifelse(cols[green]=="black","green",cols[green])
matplot(t(vc[[2]]),type="l",col=cols,lty=1,ylab="selection probability",xlab=expression(paste(lambda[v])),main="stability path columns")
lines(thr1,col="darkred",lwd=3,lty=2)
lines(thr2,col="darkorange",lwd=3,lty=2)
lines(thr3,col="darkgreen",lwd=3,lty=2)
abline(h=0.7,col="blue",lwd=2)
#abline(h=0.8,col="blue",lwd=2)
abline(v=sum(thr1<thr),lwd=2,col="darkred",lty=3)
abline(v=sum(thr2<thr),lwd=2,col="darkorange",lty=3)
abline(v=sum(thr3<thr),lwd=2,col="darkgreen",lty=3)
legend(-3, 1, c("PCER 0.1", "PCER 0.25", "PCER 0.5"),text.col = c("darkred","darkorange","darkgreen"),bty="n")
thr <- 0.7
thr1 <- ((uc[[4]]^2/((p*0.1)*p))+1)/2
red <- which(apply(uc[[2]][,1:sum(thr1<thr)],1,max)>thr)
thr2 <- ((uc[[4]]^2/((p*0.25)*p))+1)/2
orange <- which(apply(uc[[2]][,1:sum(thr2<thr)],1,max)>thr)
thr3 <- ((uc[[4]]^2/((p*0.5)*p))+1)/2
green <- which(apply(uc[[2]][,1:sum(thr3<thr)],1,max)>thr)
cols <- rep("black",nrow(uc[[2]]))
cols[red] <- "red"
cols[orange] <- ifelse(cols[orange]=="black","orange",cols[orange])
cols[green] <- ifelse(cols[green]=="black","green",cols[green])
matplot(t(uc[[2]]),type="l",col=cols,lty=1,ylab="selection probability",xlab=expression(paste(lambda[u])),main="stability path rows")
thr1 <- ((uc[[4]]^2/((p*0.1)*p))+1)/2
thr2 <- ((uc[[4]]^2/((p*0.25)*p))+1)/2
thr3 <- ((uc[[4]]^2/((p*0.5)*p))+1)/2
lines(thr1,col="darkred",lwd=3,lty=2)
lines(thr2,col="darkorange",lwd=3,lty=2)
lines(thr3,col="darkgreen",lwd=3,lty=2)
abline(h=0.7,col="blue",lwd=2)
#abline(h=0.8,col="blue",lwd=2)
abline(v=sum(thr1<thr),lwd=2,col="darkred",lty=3)
abline(v=sum(thr2<thr),lwd=2,col="darkorange",lty=3)
abline(v=sum(thr3<thr),lwd=2,col="darkgreen",lty=3)
legend(-7.5, 1, c("PCER 0.1", "PCER 0.25", "PCER 0.5"),text.col = c("darkred","darkorange","darkgreen"),bty="n")
title("Stability Paths",outer=T)





setwd("~/workspace/")
pdf("stabpath.pdf",width=16,height=9)
par(mfrow=c(1,2),omi=c(0.25, 0.25, 0.5, 0.25)) #c(bottom, left, top, right)
thr <- 0.7
thr1 <- ((vc[[4]]^2/((n*0.1)*n))+1)/2
red <- which(apply(vc[[2]][,1:sum(thr1<thr)],1,max)>thr)
thr2 <- ((vc[[4]]^2/((n*0.25)*n))+1)/2
orange <- which(apply(vc[[2]][,1:sum(thr2<thr)],1,max)>thr)
thr3 <- ((vc[[4]]^2/((n*0.5)*n))+1)/2
green <- which(apply(vc[[2]][,1:sum(thr3<thr)],1,max)>thr)
cols <- rep("black",nrow(vc[[2]]))
cols[red] <- "red"
cols[orange] <- ifelse(cols[orange]=="black","orange",cols[orange])
cols[green] <- ifelse(cols[green]=="black","green",cols[green])
matplot(t(vc[[2]]),type="l",col=cols,lty=1,ylab="selection probability",xlab=expression(paste(lambda[v])),main="stability path columns")
lines(thr1,col="darkred",lwd=3,lty=2)
lines(thr2,col="darkorange",lwd=3,lty=2)
lines(thr3,col="darkgreen",lwd=3,lty=2)
abline(h=0.7,col="blue",lwd=2)
#abline(h=0.8,col="blue",lwd=2)
abline(v=sum(thr1<thr),lwd=2,col="darkred",lty=3)
abline(v=sum(thr2<thr),lwd=2,col="darkorange",lty=3)
abline(v=sum(thr3<thr),lwd=2,col="darkgreen",lty=3)
legend(-3, 1, c("PCER 0.1", "PCER 0.25", "PCER 0.5"),text.col = c("darkred","darkorange","darkgreen"),bty="n")
thr <- 0.7
thr1 <- ((uc[[4]]^2/((p*0.1)*p))+1)/2
red <- which(apply(uc[[2]][,1:sum(thr1<thr)],1,max)>thr)
thr2 <- ((uc[[4]]^2/((p*0.25)*p))+1)/2
orange <- which(apply(uc[[2]][,1:sum(thr2<thr)],1,max)>thr)
thr3 <- ((uc[[4]]^2/((p*0.5)*p))+1)/2
green <- which(apply(uc[[2]][,1:sum(thr3<thr)],1,max)>thr)
cols <- rep("black",nrow(uc[[2]]))
cols[red] <- "red"
cols[orange] <- ifelse(cols[orange]=="black","orange",cols[orange])
cols[green] <- ifelse(cols[green]=="black","green",cols[green])
matplot(t(uc[[2]]),type="l",col=cols,lty=1,ylab="selection probability",xlab=expression(paste(lambda[u])),main="stability path rows")
thr1 <- ((uc[[4]]^2/((p*0.1)*p))+1)/2
thr2 <- ((uc[[4]]^2/((p*0.25)*p))+1)/2
thr3 <- ((uc[[4]]^2/((p*0.5)*p))+1)/2
lines(thr1,col="darkred",lwd=3,lty=2)
lines(thr2,col="darkorange",lwd=3,lty=2)
lines(thr3,col="darkgreen",lwd=3,lty=2)
abline(h=0.7,col="blue",lwd=2)
#abline(h=0.8,col="blue",lwd=2)
abline(v=sum(thr1<thr),lwd=2,col="darkred",lty=3)
abline(v=sum(thr2<thr),lwd=2,col="darkorange",lty=3)
abline(v=sum(thr3<thr),lwd=2,col="darkgreen",lty=3)
legend(-7.5, 1, c("PCER 0.1", "PCER 0.25", "PCER 0.5"),text.col = c("darkred","darkorange","darkgreen"),bty="n")
title("Stability Paths",outer=T)
dev.off()

