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
setClass('BCs4vd',
		contains = 'BiclustMethod',
		prototype = prototype(
				biclustFunction = function(X,...){s4vd(X,...)}))

BCs4vd <- function() {
	return(new('BCs4vd'))
}


#generate simulation data
dat <- datagen(meanBC=1,nG=100,nC=50,BCG=10,BCC=5,sd=1)
X <- dat[[1]]

set.seed(12345)
t1 <- system.time(res1 <- biclust(X,BCs4vd,pcer=0.05,iter=100,nbiclust=1,size=0.632,ss.thr = c(0.7,0.8),merr=0.0001,gamm=0))
set.seed(12345)
t2 <- system.time(res2 <- biclust(X,BCs4vd,pcer=0.05,iter=100,nbiclust=1,size=0.632,ss.thr = c(0.7,0.8),merr=0.0001,gamm=0,pointwise=F))
set.seed(12345)
t3 <- system.time(res3 <- biclust(X,BCs4vd,pcer=0.1,nbiclust=1,size=0.632,ss.thr = c(0.7,0.8),merr=0.001,gamm=0,pointwise=F,fullpath=TRUE))

t1
res1
t2
res2
t3
res3

stabpath(res2,1)
stabpath(res3,1)


#ependymoma data
par(mfrow=c(1,1))
X <- t(scale(t(HDsel)))

set.seed(12345)
t1 <- system.time(res1 <-  biclust(X,BCs4vd,pcer=0.1,iter=100,steps=500,nbiclust=10,size=0.632,ss.thr = c(0.6,0.7),merr=0.01,gamm=0,col.overlap=FALSE))
library(marray)
cols <- maPalette(low="blue",mid="white",high="red", k=50)
heatmapBC(X,res1,order=T,outside=T,col=cols,axes=F)

set.seed(12345)
t2 <- system.time(res2 <-  biclust(X,BCs4vd,pcer=0.3,iter=100,steps=1000,nbiclust=10,size=0.632,ss.thr = c(0.6,0.7),merr=0.01,gamm=0,col.overlap=FALSE))
heatmapBC(X,res2,number=1:res2@Number,order=T,outside=T,col=cols)

set.seed(12345)
t2 <- system.time(res3 <- s4vd(X,pcer=0.3,nbiclust=10,ss.thr = c(0.51,0.65),merr=0.05,gamm=0,col.overlap=FALSE))
heatmapBC(X,res3,order=T,outside=T,col=cols)










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






debug(s4vd)



X <- dat[[1]]
steps <- 100
size <- 0.5
pcer <- 0.2
ss.thr <- c(0.7,0.8)
gamm <- 1
iter <- 20
verbose = 1
merr = 0.001
cols.nn=TRUE
rows.nn=FALSE


startX <- X
number <- FALSE
#us <- list(NULL,NULL,FALSE)
p <- nrow(X)
n  <- ncol(X)
rowsin <- rep(TRUE,nrow(X))
colsin <- rep(TRUE,ncol(X))
stop <- FALSE
Rspath <- Cspath <- Rows <- Cols <- ul <- vl <- thrs <- list()
dl <- niter <- numeric()


SVD <- svd(X,nu=1,nv=1)
v0 <- SVD$v
u0 <- SVD$u
d0 <- SVD$d
for(i in 1:iter){
	vc <- updatev(X,u0,pcer,ss.thr,steps,size,gamm,cols.nn)
	v1 <- vc[[1]]/sqrt(sum(vc[[1]]^2)) 
	v1[is.na(v1)] <- 0
	uc <- updateu(X,v1,pcer,ss.thr,steps,size,gamm,rows.nn)
	u1 <- uc[[1]]/sqrt(sum(uc[[1]]^2)) 
	u1[is.na(u1)] <- 0
	ud <- sqrt(sum((u0-u1)^2))
	vd <- sqrt(sum((v0-v1)^2))
	if(verbose > 0) cat("iter: ",i," rows: ",length(which(u1!=0))," cols: ",length(which(v1!=0))," merr: ",min(c(ud,vd)),"\n")
	r.in <- which(u1!=0)
	c.in <- which(v1!=0)
	v0 <- v1
	u0 <- u1
	if(min(c(vd,ud)) < merr)break
}
d1 <- as.numeric(t(u0)%*%X%*%v0)	


#update u
updateu <- function(X,v0,pcer,ss.thr,steps,size,gamm,rows.nnc=FALSE,fullpath=FALSE){
	p <- nrow(X)
	err <- pcer*p
	ols <- X%*%v0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
	qs <- numeric(length(lambdas)) 
	if(rows.nnc){
		for(l in 1:length(lambdas)){
			temp <- adaLasso.nn(X,v0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thr <- ((qs[l]^2/(err*p))+1)/2
			if(thr>=ss.thr[1]&!fullpath)break
		}
	}else{
		for(l in 1:length(lambdas)){
			temp <- adaLasso(X,v0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thr <- ((qs[l]^2/(err*p))+1)/2
			if(thr>=ss.thr[1]&!fullpath)break
		}
	}
	thr <- ((qs^2/(err*p))+1)/2
	l <- sum(thr < ss.thr[1])
	thr <- thr[l]
	stable <- which(rowMeans(temp!=0)>=thr)
	if(length(stable)==0)stop <- TRUE
	uc <- rep(0,p)
	delta <- lambdas[l]/(abs(ols)^gamm)  
	uc[stable] <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))[stable]
	return(list(uc=uc,selprobpath=selprobpath,stop=stop,qs=qs,thr=thr,l=l))
}

#update u pointwise
updateu.pw <- function(X,v0,pcer,ss.thr,steps,size,gamm,rows.nnc=FALSE,l=NULL){
	p <- nrow(X)
	err <- pcer*p
	ols <- X%*%v0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
	qs <- numeric(length(lambdas)) 
	if(is.null(l)) l <- which(lambdas==quantile(lambdas,0.5,type=1))[1]
	#search for a lambda
	if(rows.nnc){
		for(g in 1:(length(lambdas))){
			temp <- adaLasso(X,v0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thr <- ((qs[l]^2/(err*p))+1)/2
			if(thr >= ss.thr[1] & thr <= ss.thr[2])break 
			if(thr < ss.thr[1]) l <- min(length(lambdas),l + ceiling(length(lambdas)/(g+1))) 
			if(thr > ss.thr[2]) l <- max(1,l - ceiling(length(lambdas)/(g+1))) 
		}
	}else{
		for(g in 1:(length(lambdas))){
			temp <- adaLasso.nn(X,v0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thr <- ((qs[l]^2/(err*p))+1)/2
			if(thr >= ss.thr[1] & thr <= ss.thr[2])break 
			if(thr < ss.thr[1]) l <- min(length(lambdas),l + ceiling(length(lambdas)/(g+1))) 
			if(thr > ss.thr[2]) l <- max(1,l - ceiling(length(lambdas)/(g+1))) 
		}
	}
	stable <- which(rowMeans(temp!=0)>=thr)
	if(length(stable)==0)stop <- TRUE
	uc <- rep(0,p)
	delta <- lambdas[l]/(abs(ols)^gamm)  
	uc[stable] <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))[stable]
	return(list(uc=uc,selprobpath=selprobpath,stop=stop,qs=qs,thr=thr,l=l))
}





system.time(res1 <- updatev(X,u0,pcer,ss.thr,steps,size,gamm,cols.nnc=T))
system.time(res2 <- updatev.pw(X,u0,pcer,ss.thr,steps,size,gamm,cols.nnc=T,l=NULL))

thr1 <- ((vc[[4]]^2/((n*pcer)*n))+1)/2
thr2 <- ((qus^2/(20*n))+1)/2
thr3 <- ((qus^2/(30*n))+1)/2
thr4 <- ((qus^2/(40*n))+1)/2
thr5 <- ((qus^2/(50*n))+1)/2
thr6 <- ((qus^2/(100*n))+1)/2
matplot(t(vc[[2]]),type="l",col="black",lty=1)
lines(thr1,col="blue",lwd=3)
lines(thr2,col="blue",lwd=3)
lines(thr3,col="blue",lwd=3)
lines(thr4,col="blue",lwd=3)
lines(thr5,col="blue",lwd=3)
lines(thr6,col="blue",lwd=3)
abline(h=ss.thr[1],col="red",lwd=3)
abline(h=0.7,col="red",lwd=3)


abline(v=sum(qus<qumax),col="red",lwd=3)
#updatev
ols <- t(X)%*%u0
lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
#lambdas <- seq(max(abs(ols)),0,length.out = p+1)
l <- 50

temp <- adaLasso(t(X),u0,lambda[l],steps,size,gamm)
qv <- mean(colSums(temp!=0))
(((vc[[4]]^2/(10*n))+1)/2)[19]

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

