# TODO: Add comment
# 
# Author: martin
###############################################################################
rm(list=ls())
dat <- datagen(noise=0.5)
X <- dat[[1]]


s4vd <- function(
		steps <- 100
		size <- 0.5
		pcer <- 0.2
		ss.thr <- c(0.7,0.8)
		gamm <- 1
		iter <- 20
		verbose = 1
		merr = 0.05
		cols.nn=TRUE
		rows.nn=FALSE
		row.overlap=TRUE
		col.overlap=TRUE
){
	startX <- X
	number <- FALSE
	us <- list(NULL,NULL,FALSE)
	p <- nrow(X)
	n <- ncol(X)
	rowsin <- rep(TRUE,p)
	colsin <- rep(TRUE,n)
	stop <- FALSE
	info <- Rows <- Cols <- list()
	
	for(k in 1:nbiclust){
		cat("Bicluster",k,"\n")
		rows <- rep(FALSE,nrow(startX))
		cols <- rep(FALSE,ncol(startX))
		if(is.null(nrow(X))|is.null(ncol(X))){
			break
		}
		if(nrow(X)==0|ncol(X)==0){
			break
		}
		SVD <- svd(X,nu=1,nv=1)
		v0 <- SVD$v
		u0 <- SVD$u
		d0 <- SVD$d
		if((length(u0)*size)<2|(length(v0)*size)<=2){
			cat("submatrix to small for resampling","\n")
			number <- k-1
			stop <- TRUE
			break
		}
		for(i in 1:iter){
			vc <- updatev(X,u0,pcer,ss.thr,steps,size,gamm,cols.nn,fullpath=T)
			v1 <- vc[[1]]/sqrt(sum(vc[[1]]^2)) 
			v1[is.na(v1)] <- 0
			uc <- updateu(X,v1,pcer,ss.thr,steps,size,gamm,rows.nn,fullpath=T)
			u1 <- uc[[1]]/sqrt(sum(uc[[1]]^2)) 
			u1[is.na(u1)] <- 0
			ud <- sqrt(sum((u0-u1)^2))
			vd <- sqrt(sum((v0-v1)^2))
			cat("iter: ",i," rows: ",length(which(u1!=0))," cols: ",length(which(v1!=0))," merr: ",min(c(ud,vd)),"\n")
			r.in <- which(u1!=0)
			c.in <- which(v1!=0)
			v0 <- v1
			u0 <- u1
			if(min(c(vd,ud)) < merr)break
		}	
		rows[rowsin] <- u0!=0
		cols[colsin] <- v0!=0
		Rows[[k]] <- rows
		Cols[[k]] <- cols
		if(!row.overlap){
			rowsin[rows] <- FALSE
			X <- startX[rowsin,colsin]
			} 
		if(!col.overlap){
			colsin[cols] <- FALSE
			X <- startX[rowsin,colsin]
			} 
		if(row.overlap&col.overlap){
			d0 <- as.numeric(t(u0)%*%X%*%v0)
			X <- round(X - d0*u0%*%t(v0))
			}
		if(stop){
			number <- k-1
			cat("no further stable bicluster, increase PCER!","\n")
			break
		}
		if(i==iter){
			number <- k-1
			cat("bicluster not coverged, increase merr or the number of iterations!","\n")
			break
		}
	}
	if(!number) number <- k
	niter[k] <- e
	Rows[[k]] <- rows
	Cols[[k]] <- cols
	params <- list(steps = steps, pcer=pcer, iter=iter, ss.thr=ss.thr, size=size, gamm=gamm, r.overlap=r.overlap, c.overlap=c.overlap,
			r.negcorr=r.negcorr, c.negcorr=c.negcorr, nbiclust=nbiclust, merr=merr)
	RowxNumber=t(matrix(unlist(Rows),byrow=T,ncol=length(Rows[[1]])))
	NumberxCol=matrix(unlist(Cols),byrow=T,ncol=length(Cols[[1]]))
	Number <- number
	info <- list(ul=ul,vl=vl,dl=dl,Rspath=Rspath,Cspath=Cspath,niter=niter,thrs=thrs)
	return(BiclustResult(params,RowxNumber,NumberxCol,Number,info))
}


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
	return(list(uc,selprobpath,stop,qs,thr,l))
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
	return(list(uc,selprobpath,stop,qs,thr,l))
}

#update v
updatev <- function(X,u0,pcer,ss.thr,steps,size,gamm,cols.nnc=FALSE,fullpath=FALSE){
	n <- ncol(X)
	err <- pcer*n
	ols <- t(X)%*%u0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
	qs <- numeric(length(lambdas)) 
	if(cols.nnc){
		for(l in 1:length(lambdas)){
			temp <- adaLasso.nn(t(X),u0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thr <- ((qs[l]^2/(err*n))+1)/2
			if(thr>=ss.thr[1]&!fullpath)break
		}
	}else{
		for(l in 1:length(lambdas)){
			temp <- adaLasso(t(X),u0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thr <- ((qs[l]^2/(err*n))+1)/2
			if(thr>=ss.thr[1]&!fullpath)break
		}
	}
	thr <- ((qs^2/(err*n))+1)/2
	l <- sum(thr < ss.thr[1])
	thr <- thr[l]
	stable <- which(rowMeans(temp!=0)>=thr)
	if(length(stable)==0)stop <- TRUE
	vc <- rep(0,n)
	delta <- lambdas[l]/(abs(ols)^gamm)  
	vc[stable] <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))[stable]
	return(list(vc,selprobpath,stop,qs,thr))
}

updatev.pw <- function(X,u0,pcer,ss.thr,steps,size,gamm,cols.nnc=FALSE,l=NULL){
	n <- ncol(X)
	err <- pcer*n
	ols <- t(X)%*%u0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
	qs <- numeric(length(lambdas)) 
	if(is.null(l)) l <- which(lambdas==quantile(lambdas,0.5,type=1))[1]
	#search for a lambda
	if(cols.nnc){
		for(g in 1:(length(lambdas))){
			temp <- adaLasso(t(X),u0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thr <- ((qs[l]^2/(err*n))+1)/2
			if(thr >= ss.thr[1] & thr <= ss.thr[2])break 
			if(thr < ss.thr[1]) l <- min(length(lambdas),l + ceiling(length(lambdas)/(g+1))) 
			if(thr > ss.thr[2]) l <- max(1,l - ceiling(length(lambdas)/(g+1))) 
		}
	}else{
		for(g in 1:(length(lambdas))){
			temp <- adaLasso.nn(t(X),u0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thr <- ((qs[l]^2/(err*n))+1)/2
			if(thr >= ss.thr[1] & thr <= ss.thr[2])break 
			if(thr < ss.thr[1]) l <- min(length(lambdas),l + ceiling(length(lambdas)/(g+1))) 
			if(thr > ss.thr[2]) l <- max(1,l - ceiling(length(lambdas)/(g+1))) 
		}
	}
	stable <- which(rowMeans(temp!=0)>=thr)
	if(length(stable)==0)stop <- TRUE
	vc <- rep(0,n)
	delta <- lambdas[l]/(abs(ols)^gamm)  
	vc[stable] <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))[stable]
	return(list(vc,selprobpath,stop,qs,thr,l))
}


adaLasso <- function(X,b,lambda,steps,size,gamm=0){
	subsets <- sapply(1:steps,function(x){sample(1:length(b),length(b)*size)})
	res <- sapply(1:steps,adaLassoSteps,subsets,X,b,lambda,gamm)
	return(res)
}
adaLassoSteps <- function(index,subsets,X,b,lambda,gamm){
	ols <- X[,subsets[,index]]%*%b[subsets[,index]]
	delta <- lambda/(abs(ols)^gamm)                        
	ols <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	ols[is.na(ols)] <- 0
	return(ols)
}
adaLasso.nn <- function(X,b,lambda,steps,size,gamm=0){
	subsets <- sapply(1:steps,function(x){sample(1:length(b),length(b)*size)})
	res <- sapply(1:steps,adaLassoSteps.nn,subsets,X,b,lambda,gamm)
	return(res)
}
adaLassoSteps.nn <- function(index,subsets,X,b,lambda,gamm){
	ols <- X[,subsets[,index]]%*%b[subsets[,index]]
	delta <- lambda/(abs(ols)^gamm)                        
	ols <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	ols[is.na(ols)] <- 0
	mostof <- sign(sum(sign(ols)))
	if(mostof==0) mostof <- 1
	ols[which(sign(ols) != mostof) ] <- 0
	ols[is.na(ols)] <- 0
	return(ols)
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


par(mfrow=c(1,2))
matplot(t(vc[[2]]),type="l",col="black",lty=1)
thr1 <- ((vc[[4]]^2/((n*pcer)*n))+1)/2
thr2 <- ((vc[[4]]^2/((n*0.3)*n))+1)/2
thr3 <- ((vc[[4]]^2/((n*0.5)*n))+1)/2
lines(thr1,col="blue",lwd=3)
lines(thr2,col="blue",lwd=3)
lines(thr3,col="blue",lwd=3)
matplot(t(uc[[2]]),type="l",col="black",lty=1)
thr1 <- ((uc[[4]]^2/((n*pcer)*p))+1)/2
thr2 <- ((uc[[4]]^2/((n*0.3)*p))+1)/2
thr3 <- ((uc[[4]]^2/((n*0.5)*p))+1)/2
lines(thr1,col="blue",lwd=3)
lines(thr2,col="blue",lwd=3)
lines(thr3,col="blue",lwd=3)



