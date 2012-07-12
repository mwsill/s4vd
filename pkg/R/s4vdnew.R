rm(list=ls())
library(parallel)
library(s4vd)
library(Matrix)

data(lung) 
X <- lung


plot(density(lambdas))

regPath <- function(X,v,gamm,lambdas){
  ols <- as.numeric(X%*%v)
  lambdaMat <- t(matrix(rep(lambdas,length(ols)),ncol=length(ols)))
  olsMat    <- matrix(rep(ols,length(lambdas)),ncol=length(lambdas))
  lambdaMat <- lambdaMat/(abs(olsMat)^gamm)
  regPath <- sign(olsMat) * (abs(olsMat) >= lambdaMat) * (abs(olsMat) - lambdaMat) 
  regPath != 0
}

stabPath <- function(X,v,gamm=gamm,size=size,steps=steps,cores=mc.cores){
	n <- ncol(X)
	p <- nrow(X)
	subsets <- sapply(1:steps,function(x){sample(1:n,n*size)})
	ols <- as.numeric(X %*% v)
	lambda.min.ratio = ifelse(n < p, 0.01,1e-04)
	#length of the path is limited to 100 lambdas
	lambdas <- exp(seq(log(max(abs(ols))),log(max(abs(ols))*lambda.min.ratio),length.out=100))
	system.time(
	res <- mclapply(1:steps,mc.cores=cores,function(index)
	       regPath(X[,subsets[,index]],v[subsets[,index]],0,lambdas)))
	#merging
	stabpath <- Reduce("+",res)/steps
	qmat <- sapply(res,function(x)colSums(x))
	rm(res)
	gc()
	qs <- rowMeans(qmat)
	return(list(stabpath=stabpath,qs=qs,lambdas=lambdas))
}

stabilitySelection <- function(X,v,pcer,pi_thr,steps,size,mc.cores,gamm){
	spath <- stabPath(X,v,gamm=gamm,size=size,steps=steps,cores=mc.cores)
	p <- nrow(X)
	qthr <- sqrt(pi_thr*pcer*p*p)
	lpos <- which(spath[[2]]>qthr)[1]
	stable <- spath[[1]][,lpos]>=pi_thr
	ols <- X%*%v
	lambda <- spath[[3]][lpos]
	lambda <- lambda/(abs(ols)^gamm)
	#soft thresholding			
	hardOLS <- softOLS <- sign(ols) * (abs(ols) >= lambda) * (abs(ols) - lambda)
	#softOLS <- sign(ols) * (abs(ols) >= lambda) * (abs(ols) - lambda)
	#hardOLS <- ols
	#hard thresholding of all variables not included in the stable set 
	hardOLS[!stable] <- 0	
	return(list(softOLS=softOLS,hardOLS=hardOLS,spath=spath,lambda=lambda,lpos=lpos))
}

s4vdLayer <- function(X,pceru=.1,pcerv=NULL,pi.thr=0.8,b=100,max.iter=20,size=0.632
                      ,mc.cores=getOption("mc.cores", 2L),gamm=0,merr=10^(-3)){
	  SVD <- svd(X, nu = 1, nv = 1)
    v <- SVD$v
    u <- SVD$u
    for(i in 1:max.iter){
		  uweights <- rep(gamm,length(u))
		  vweights <- rep(gamm,length(v))
		  resu <- resv <- list()
		  #update u
		  if(!is.null(pceru)){
		  	resu <-stabilitySelection(X,v,pceru,pi_thr,
			  		steps,size,mc.cores,gamm=uweights)
			u1 <- as.numeric(resu$softOLS)
			#uweights <- 1-resu$spath[[1]][,res$lpos]
		  }else{
			  u1 <- X%*%v
		  }
		  u1 <- u1/sqrt(sum(u1^2))
		  u1[is.na(u1)] <- 0
		  #update v
		  if(!is.null(pcerv)){
			resv <- stabilitySelection(t(X),u1,pcerv,pi_thr,
					steps,size,mc.cores,gamm=vweights)
			v1 <- as.numeric(resv$softOLS)
			#vweights <- 1-resv$spath[[1]][,resv$lpos]
		  }else{
			  v1 <- t(X)%*%u1
		  }
		  v1 <- v1/sqrt(sum(v1^2))
		  v1[is.na(v1)] <- 0
		  # check convergence	
		  ud <- sqrt(sum((u - u1)^2))
		  vd <- sqrt(sum((v - v1)^2))
		  #cat("ud",ud,"vd",vd,"\n")
		  cat(".")
      v <- v1
  		u <- u1
  		if(min(c(ud,vd))<merr) break
  	}
  	if(!is.null(pcerv)){ 
      v <- resv$hardOLS
  		v1 <- v1/sqrt(sum(v1^2))
  		v1[is.na(v1)] <- 0
  	} 
  	if(!is.null(pceru)){ 
      u <- resu$hardOLS
  		u <- u/sqrt(sum(u^2))
  		u[is.na(u)] <- 0
  	}
  	cat("\n") 	
  	d <- as.numeric(t(u) %*% X %*% v)	#as.numeric(t(u0) %*% X %*% v0)
  	return(list(u=u,v=v,d=d,selectionU=resu,selectionV=resv))	
}
res <- s4vdLayer(X,pceru=.5,pcerv=NULL,pi_thr=.8,steps=100,max.iter=20,size=.632,mc.cores=2,gamm=0,merr=0.01)


#s4vdBiclustering <- function(X,fweru=.5,fwerv=.5,pi_thr=.8,steps=100,max.iter=20,size=.632
#		,mc.cores=8,gamm=0,row.overlap=TRUE,col.overlap=TRUE,ncluster=5){
#		cluster <- list()
#		or(i in 1:ncluster){
#		res <- s4vdLayer(X,fweru,fwerv,pi_thr,steps,max.iter,size,mc.cores,gamm)
#		Xnew <- X
#		}	
#}

s4vd <- function(X,k=3,pceru=.1,pcerv=.1,pi.thr=.8,b=100
,max.iter=20,size=.632,mc.cores=8,gamm=0,merr=10^(-3)){
	MYCALL <- match.call()
	u <- matrix(nrow=nrow(X),ncol=k)
	v <- matrix(nrow=ncol(X),ncol=k)
	d <- numeric()
	for(i in 1:k){
  	cat("Layer",i)
  	temp <- s4vdLayer(X,pceru=pceru,pcerv=pcerv,pi.thr,b,max.iter,size,mc.cores,gamm,merr)
  	u[,i] <- temp[[1]] 
  	v[,i] <- temp[[2]]
  	d[i] <-  temp[[3]]
  	X <- X - d[i]*u[,i]%*%t(v[,i])
	}	
	return(list(u=u,v=v,d=d))
}




#parallel implementation via multicore
system.time(
res <- s4vdPCA(X,K=3,pcer=.1,pi_thr=.8,steps=100
	      ,size=.632,max.iter=20,mc.cores=1))

#parallel implementation via multicore 
system.time(
res <- s4vdPCA(X,K=3,pcer=.1,pi_thr=.8,steps=100
	      ,size=.632,max.iter=20,mc.cores=8))

#plot the stability path for PC1
stabpath(res,1)
#biplot 2D and 3D
biplot(res,d3=FALSE)

X <- lung + matrix(rnorm(nrow(X)*ncol(X),0,1.5),nrow=nrow(X),ncol=ncol(X))

res <-  svd(X) 
#res1 <- res
#biplot
#par(mfrow=c(1,1),omi=c(0,0,0,2))
alpha <- 0.4

par(mfrow=c(2,2))
ptx <- res$v[,1]*(res$d[1]^alpha)
pty <- res$v[,2]*(res$d[1]^alpha)
ptz <- res$v[,3]*(res$d[1]^alpha)
arrowsx <- res$u[,1]*(res$d[1]^(1-alpha))
arrowsy <- res$u[,2]*(res$d[1]^(1-alpha))
arrowsz <- res$u[,3]*(res$d[1]^(1-alpha))
obj.cols <-ifelse(colnames(X)=="Carcinoid","orange","darkgreen")
obj.cols[colnames(X)=="Colon"] <- "blue"
obj.cols[colnames(X)=="SmallCell"] <- "darkred"



plot(ptx,pty,type="n",xlim=c(-5,5),ylim=c(-5,5),axes=T,ylab="PC2",xlab="PC1",main="Biplot sparse PCA",lty=2,col=obj.cols)
arrows(0, 0, x1 = arrowsx, y1 =arrowsy ,length=0.05,col = "pink")
points(ptx,pty,col=obj.cols)

plot(ptx,ptz,type="n",xlim=c(-5,5),ylim=c(-5,5),axes=T,ylab="PC3",xlab="PC1",main="Biplot sparse PCA",lty=2,col=obj.cols)
arrows(0, 0, x1 = arrowsx, y1 =arrowsz ,length=0.05,col = "pink")
points(ptx,ptz,col=obj.cols)

plot(pty,ptz,type="n",xlim=c(-5,5),ylim=c(-5,5),axes=T,ylab="PC2",xlab="PC3",main="Biplot sparse PCA",lty=2,col=obj.cols)
arrows(0, 0, x1 = arrowsy, y1 =arrowsz ,length=0.05,col = "pink")
points(pty,ptz,col=obj.cols)

#rotating 3d biplot
library(rgl)

plot3d(ptx,pty,ptz,xlim=c(min(c(min(ptx),min(arrowsx))),max(c(max(ptx),max(arrowsx)))),
             ylim=c(min(c(min(pty),min(arrowsy))),max(c(max(pty),max(arrowsy)))),
             zlim=c(min(c(min(ptz),min(arrowsz))),max(c(max(ptz),max(arrowsz)))),
             xlab='',
             ylab='',
             zlab='',
             type='n',
             axes=FALSE,
             box=box,
             aspect=TRUE,
             top=TRUE)

spheres3d(ptx,pty,ptz,
                  col=obj.cols,
                  radius=.25,
                  alpha=.75)

vars<- cbind(arrowsx,arrowsy,arrowsz)
apply(vars,1, function(x)segments3d(rbind(c(0,0,0),x),col="red",alpha=.25))

movie3d( spin3d(), duration=15 ,dir="/windows/workspace/useR2012/vortrag2.0/ani/PCA",clean=FALSE)


PC1 <- s4vdLayer(X,pceru=.05,pcerv=NULL,pi_thr=.8,steps=100,max.iter=20,size=.632,mc.cores=8,gamm=0,merr=10^(-2))


cols <- rep("gray",nrow(X))
cols[which(PC1$selectionU$hardOLS!=0)] <- "red"

pdf("stabpath.pdf")
matplot(t(PC1$selectionU$spath[[1]]),xaxt="n",type="l",col=cols,lty=1,ylab=expression(paste(hat(Pi)))
        ,xlab=expression(paste(lambda)),main="Stability Path",cex.lab=1,cex.axis=1)
abline(h=0.8,col="darkred",lwd=2)
abline(v=PC1$selectionU$lpos,col="darkred",lwd=2)
dev.off()



#res <- s4vdPCA(X,K=3,pcer=.05,mc.cores=8,merr=10^(-2))
SVD <- svd(X)
n <- ncol(X)
p <- nrow(X)
v <- SVD$v[,1]
ols <- as.numeric(X %*% v)
lambda.min.ratio = ifelse(n < p, 0.01,1e-04)
#length of the path is limited to 100 lambdas
lambdas <- exp(seq(log(max(abs(ols))),log(max(abs(ols))*lambda.min.ratio),length.out=100))
gamm <- 0
ols <- as.numeric(X%*%v)
  lambdaMat <- t(matrix(rep(lambdas,length(ols)),ncol=length(ols)))
  olsMat    <- matrix(rep(ols,length(lambdas)),ncol=length(lambdas))
  lambdaMat <- lambdaMat/(abs(olsMat)^gamm)
  regPath <- sign(olsMat) * (abs(olsMat) >= lambdaMat) * (abs(olsMat) - lambdaMat)

pdf("regpath.pdf")
matplot(t(regPath),xaxt="n",yaxt="n",type="l",col="black",lty=1,ylab=expression(paste(beta[i]))
        ,xlab=expression(paste(lambda)),main="Penalization Path",cex.lab=1,cex.axis=1)
dev.off()








alpha <- 0.4

par(mfrow=c(2,2))
ptx <- res$v[,1]*(res$d[1]^alpha)
pty <- res$v[,2]*(res$d[1]^alpha)
ptz <- res$v[,3]*(res$d[1]^alpha)
arrowsx <- res$u[,1]*(res$d[1]^(1-alpha))
arrowsy <- res$u[,2]*(res$d[1]^(1-alpha))
arrowsz <- res$u[,3]*(res$d[1]^(1-alpha))
obj.cols <-ifelse(colnames(X)=="Carcinoid","orange","darkgreen")
obj.cols[colnames(X)=="Colon"] <- "blue"
obj.cols[colnames(X)=="SmallCell"] <- "darkred"



plot(ptx,pty,type="n",xlim=c(-5,5),ylim=c(-5,5),axes=T,ylab="PC2",xlab="PC1",main="Biplot sparse PCA",lty=2,col=obj.cols)
arrows(0, 0, x1 = arrowsx, y1 =arrowsy ,length=0.05,col = "pink")
points(ptx,pty,col=obj.cols)

plot(ptx,ptz,type="n",xlim=c(-5,5),ylim=c(-5,5),axes=T,ylab="PC3",xlab="PC1",main="Biplot sparse PCA",lty=2,col=obj.cols)
arrows(0, 0, x1 = arrowsx, y1 =arrowsz ,length=0.05,col = "pink")
points(ptx,ptz,col=obj.cols)

plot(pty,ptz,type="n",xlim=c(-5,5),ylim=c(-5,5),axes=T,ylab="PC2",xlab="PC3",main="Biplot sparse PCA",lty=2,col=obj.cols)
arrows(0, 0, x1 = arrowsy, y1 =arrowsz ,length=0.05,col = "pink")
points(pty,ptz,col=obj.cols)

#rotating 3d biplot
library(rgl)

plot3d(ptx,pty,ptz,xlim=c(min(c(min(ptx),min(arrowsx))),max(c(max(ptx),max(arrowsx)))),
             ylim=c(min(c(min(pty),min(arrowsy))),max(c(max(pty),max(arrowsy)))),
             zlim=c(min(c(min(ptz),min(arrowsz))),max(c(max(ptz),max(arrowsz)))),
             xlab='',
             ylab='',
             zlab='',
             type='n',
             axes=FALSE,
             box=box,
             aspect=TRUE,
             top=TRUE)

spheres3d(ptx,pty,ptz,
                  col=obj.cols,
                  radius=.25,
                  alpha=.75)

vars<- cbind(arrowsx,arrowsy,arrowsz)
apply(vars,1, function(x)segments3d(rbind(c(0,0,0),x),col="red",alpha=.25))


movie3d( spin3d(), duration=15 ,dir="/windows/workspace/useR2012/vortrag2.0/ani/sPCA",clean=FALSE)


axes3d(c("PC1","PC2","PC3"))
        # simple axis
        title3d(xlab=xlab,
                ylab=ylab,
                zlab=zlab, ...)








plot(ptx,pty,type="p",xlim=c(min(c(min(ptx),min(arrowsx)))
,max(c(max(ptx),max(arrowsx)))),ylim=c(min(c(min(pty),min(arrowsy)))
,max(c(max(pty),max(arrowsy)))),axes=T,ylab="PC2",xlab="PC1",main="Biplot sparse PCA",lty=2)
arrows(0, 0, x1 = arrowsx, y1 =arrowsy ,length=0.05,col = "pink")
#arrows(0, 0, x1 = pc1$v[top10]*(pc1$d^(1-alpha)), y1 = pc2$v[top10]*(pc2$d^(1-alpha)),length=0.05,col = "red")



image(res$d*res$u%*%t(res$v))


par(mfrow=(1,2))
image(X)
image(res$u*res$d%*%res$v)

     	MYCALL <- match.call()
    	startX <- X
    	p.ini <- nrow(X)
    	n.ini <- ncol(X)
    	rowsin <- rep(TRUE, p.ini)
    	colsin <- rep(TRUE, n.ini)
    	stop <- FALSE
    	start <- TRUE
    	info <- Rows <- Cols <- vc <- uc <- list()
    	for (k in 1:nbiclust) {
        cat("Bicluster", k)
        rows <- rep(FALSE, nrow(startX))
        cols <- rep(FALSE, ncol(startX))
        if (is.null(nrow(X)) | is.null(ncol(X))) {
            number <- k - 1
            stop <- TRUE
            break
        }
        if (nrow(X) == 0 | ncol(X) == 0) {
            number <- k - 1
            stop <- TRUE
            break
        }
        SVD <- svd(X, nu = 1, nv = 1)
        v0 <- SVD$v
        u0 <- SVD$u
        d0 <- SVD$d
        vc <- uc <- list()
        if ((length(u0) * size) <= 2 | (length(v0) * size) <= 
            2) {
            cat("submatrix to small for resampling", "\n")
            number <- k - 1
            stop <- TRUE
            break
        }
        if (pointwise) {
            for (i in 1:iter) {
                if (i > start.iter) 
                  start <- FALSE
                uc <- updateu.pw(X, v0, pceru, p.ini, ss.thr, 
                  steps, size, gamm, rows.nc, uc$l)
                u1 <- uc[[1]]/sqrt(sum(uc[[1]]^2))
                u1[is.na(u1)] <- 0
                vc <- updatev.pw(X, u1, pcerv, n.ini, ss.thr, 
                  steps, size, gamm, cols.nc, vc$l)
                v1 <- vc[[1]]/sqrt(sum(vc[[1]]^2))
                v1[is.na(v1)] <- 0
                if (uc[[3]] & i > start.iter) {
                  cat("rows not stable")
                  stop <- TRUE
                  break
                }
                if (vc[[3]] & i > start.iter) {
                  cat("columns not stable")
                  stop <- TRUE
                  break
                }
                ud <- sqrt(sum((u0 - u1)^2))
                vd <- sqrt(sum((v0 - v1)^2))
                cat(".")
                v0 <- v1
                u0 <- u1
                if (min(c(vd, ud)) < merr & i > start.iter) 
                  break
            }
        }
}





#u1[is.na(u1)] <- 0
fwer <- 0.05
pi_thr <- .8

matplot(t(spath[[1]]),type="l",col="black",lty=1)
	abline(h=pi_thr,col="red")
	abline(v=lpos,col="red")

matplot(t(res$spathSstabpath),type="l")


           
matplot(t(regPath(t(X),u,0)),type="l")

#draw subsets
subsets <- sapply(1:steps,function(x){sample(1:p,p*size)})

uc <- rep(0, p)
delta <- lambdas[ls]/(abs(ols)^gamm)
uc <- (sign(ols) * (abs(ols) >= delta) * (abs(ols) - delta))



s4vd
function (X, steps = 100, pcerv = 0.1, pceru = 0.1, ss.thr = c(0.6, 
    0.65), size = 0.5, gamm = 0, iter = 100, nbiclust = 10, merr = 10^(-3), 
    cols.nc = TRUE, rows.nc = TRUE, row.overlap = TRUE, col.overlap = TRUE, 
    row.min = 1, col.min = 1, pointwise = TRUE, start.iter = 3, 
  
    MYCALL <- match.call()
    startX <- X
    p.ini <- nrow(X)
    n.ini <- ncol(X)
    rowsin <- rep(TRUE, p.ini)
    colsin <- rep(TRUE, n.ini)
    stop <- FALSE
    start <- TRUE
    info <- Rows <- Cols <- vc <- uc <- list()
    for (k in 1:nbiclust) {
        gc()
        cat("Bicluster", k)
        rows <- rep(FALSE, nrow(startX))
        cols <- rep(FALSE, ncol(startX))
        if (is.null(nrow(X)) | is.null(ncol(X))) {
            number <- k - 1
            stop <- TRUE
            break
        }
        if (nrow(X) == 0 | ncol(X) == 0) {
            number <- k - 1
            stop <- TRUE
            break
        }
        SVD <- svd(X, nu = 1, nv = 1)
        v0 <- SVD$v
        u0 <- SVD$u
        d0 <- SVD$d
        vc <- uc <- list()
        if ((length(u0) * size) <= 2 | (length(v0) * size) <= 
            2) {
            cat("submatrix to small for resampling", "\n")
            number <- k - 1
            stop <- TRUE
            break
        }
        if (pointwise) {
            for (i in 1:iter) {
                if (i > start.iter) 
                  start <- FALSE
                uc <- updateu.pw(X, v0, pceru, p.ini, ss.thr, 
                  steps, size, gamm, rows.nc, uc$l)
                u1 <- uc[[1]]/sqrt(sum(uc[[1]]^2))
                u1[is.na(u1)] <- 0
                vc <- updatev.pw(X, u1, pcerv, n.ini, ss.thr, 
                  steps, size, gamm, cols.nc, vc$l)
                v1 <- vc[[1]]/sqrt(sum(vc[[1]]^2))
                v1[is.na(v1)] <- 0
                if (uc[[3]] & i > start.iter) {
                  cat("rows not stable")
                  stop <- TRUE
                  break
                }
                if (vc[[3]] & i > start.iter) {
                  cat("columns not stable")
                  stop <- TRUE
                  break
                }
                ud <- sqrt(sum((u0 - u1)^2))
                vd <- sqrt(sum((v0 - v1)^2))
                cat(".")
                v0 <- v1
                u0 <- u1
                if (min(c(vd, ud)) < merr & i > start.iter) 
                  break
            }
        }
        else {
            for (i in 1:iter) {
                if (i > start.iter) 
                  start <- FALSE
                uc <- updateu(X, v0, pceru, p.ini, ss.thr, steps, 
                  size, gamm, rows.nc, savepath)
                u1 <- uc[[1]]/sqrt(sum(uc[[1]]^2))
                u1[is.na(u1)] <- 0
                vc <- updatev(X, u1, pcerv, n.ini, ss.thr, steps, 
                  size, gamm, cols.nc, savepath)
                v1 <- vc[[1]]/sqrt(sum(vc[[1]]^2))
                v1[is.na(v1)] <- 0
                if (uc[[3]] & i > start.iter) {
                  cat("rows not stable")
                  stop <- TRUE
                  break
                }
                if (vc[[3]] & i > start.iter) {
                  cat("columns not stable")
                  stop <- TRUE
                  break
                }
                ud <- sqrt(sum((u0 - u1)^2))
                vd <- sqrt(sum((v0 - v1)^2))
                cat(".")
                v0 <- v1
                u0 <- u1
                if (min(c(vd, ud)) < merr & i > start.iter) 
                  break
            }
        }
        stableu <- uc$sp >= uc$thr
        stablev <- vc$sp >= vc$thr
        d0 <- as.numeric(t(u0) %*% X %*% v0)
        u0[!stableu] <- 0
        v0[!stablev] <- 0
        rows[rowsin] <- u0 != 0
        cols[colsin] <- v0 != 0
        Rows[[k]] <- rows
        Cols[[k]] <- cols
        if (stop) {
            number <- k - 1
            break
        }
        if (i == iter) {
            number <- k - 1
            stop <- TRUE
            cat("Fail to converge! Increase the number of iterations !", 
                "\n")
            gc()
            break
        }
        if (!row.overlap) {
            rowsin[rows] <- FALSE
            X <- startX[rowsin, colsin]
            info[[k]] <- list(vc, uc, layer = list(u0, v0, d0))
        }
        if (!col.overlap) {
            colsin[cols] <- FALSE
            X <- startX[rowsin, colsin]
            info[[k]] <- list(vc, uc, layer = list(u0, v0, d0))
        }
        if (row.overlap & col.overlap) {
            temp <- svd(X[rows, cols])
            X[rows, cols] <- X[rows, cols] - (temp$d[1] * temp$u[, 
                1] %*% t(temp$v[, 1]))
            info[[k]] <- list(vc, uc, layer = list(u0, v0, d0))
        }
        cat("\n")
    }
    if (!stop) 
        number <- k
    params <- list(steps = steps, pcerv = pcerv, pceru = pceru, 
        iter = iter, ss.thr = ss.thr, size = size, gamm = gamm, 
        row.overlap = row.overlap, col.overlap = col.overlap, 
        rows.nc = rows.nc, cols.nc = cols.nc, nbiclust = nbiclust, 
        merr = merr, row.min = row.min, col.min = col.min, pointwise = pointwise, 
        start.iter = start.iter, savepath = savepath, Call = MYCALL)
    RowxNumber = t(matrix(unlist(Rows), byrow = T, ncol = length(Rows[[1]])))
    NumberxCol = matrix(unlist(Cols), byrow = T, ncol = length(Cols[[1]]))
    if (number) 
        RowxNumber <- matrix(RowxNumber[, 1:number], ncol = number)
    if (number) 
        NumberxCol <- matrix(NumberxCol[1:number, ], nrow = number)
    Number <- number
    info[[k + 1]] <- params
    return(BiclustResult(params, RowxNumber, NumberxCol, Number, 
        info))
}

