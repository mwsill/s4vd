s4vd <- function(
		X        		    # the p times n data matrix (p>>n)
		,steps=100          # number of subsamples used to calculate the relative selection frequencies for left singular vector coefficients that correspond to the rows of X
		,pcer=0.05          # expected number of falsly selected left singular vector coefficients 	
		,iter=50	        # maximal number of iterations to fit a single bicluster	
		,ss.thr=0.6         # cutoff threshold (relative selection frequency) for the stability selection
		,size=0.5           # size of the subsamples for the stability selection
		,weak=0.6           # weakness parameter for the randomised lasso
		,r.overlap=TRUE     # allow bicluster to be column overlapping
		,c.overlap=TRUE     # allow bicluster to be row overlapping		
		,r.negcorr=TRUE     # allow for negative correlation of rows (genes) over columns (samples)  
		,c.negcorr=FALSE    # allow for negative correlation of columns (samples) over rows (genes)
		,nbiclust=20        # maximal number of biclusters
		,merr=0.1           # convergence parameter 
		,mc.cores=1         # number of cores for parallelization, computes the complete path
){
	startX <- X
	number <- FALSE
	p <- nrow(X)
	n  <- ncol(X)
	rowsin <- rep(TRUE,nrow(X))
	colsin <- rep(TRUE,ncol(X))
	stop <- FALSE
	Rspath <- Cspath <- Rows <- Cols <- ul <- vl <- list()
	dl <- niter <- numeric()
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
		Xsvd <- svd(X)
		u0 <- u.ini <- Xsvd$u[,1]
		v0 <- v.ini <- Xsvd$v[,1]
		d0 <- d.ini <- Xsvd$d[1] 
		if((length(u0)*size)<2|(length(v0)*size)<=2){
			cat("submatrix to small for resampling","\n")
			number <- k-1
			stop <- TRUE
			break
		}
		for(e in 1:iter){
			#cat(".")
			vs <- updatev(X, u0, steps, c.err= pcer*ncol(X), ss.thr,size, weak, c.negcorr, mc.cores)
			v1 <- vs[[1]]/sqrt(sum(vs[[1]]^2)) 
			v1[is.na(v1)] <- 0
			cat("cols: ",length(which(v1!=0)),"\n")
			if(vs[[3]]){
				stop <- TRUE
				break
			}
			us <- updateu(X, v1, steps, r.err= pcer*nrow(X), ss.thr,size, weak, r.negcorr, mc.cores)
			u1 <- us[[1]]/sqrt(sum(us[[1]]^2))
			u1[is.na(u1)] <- 0
			cat("rows: ",length(which(u1!=0)),"\n")
			if(us[[3]]){
				stop <- TRUE
				break
			}
			ud <- sqrt(sum((u0-u1)^2))
			vd <- sqrt(sum((v0-v1)^2))
			vm <- mean(abs(u0-u1))
			um <- mean(abs(v0-v1))
			cat("iter: ",e," merr: ",mean(c(ud,vd)),"\n")
			cat("iter: ",e," merr: ",mean(c(um,vm)),"\n")
			r.in <- which(u1!=0)
			c.in <- which(v1!=0)
			d1 <- as.numeric(t(u0)%*%X%*%v0)
			v0 <- v1
			u0 <- u1
			d0 <- d1
			#cat("iter: ",e," merr: ",min(vd,ud),"\n")
			if(mean(c(vd,ud)) < merr){
				cat("\n")
				niter[k] <- e
				Rspath[[k]] <- us[[2]]
				Cspath[[k]] <- vs[[2]]
				rows[rowsin] <- u0!=0
				cols[colsin] <- v0!=0
				Rows[[k]] <- rows
				Cols[[k]] <- cols
				ul[[k]] <- u0
				vl[[k]] <- v0
				dl[k] <- d0
				if(!r.overlap){
					rowsin[rows] <- FALSE
					X <- startX[rowsin,colsin]
				} 
				if(!c.overlap){
					colsin[cols] <- FALSE
					X <- startX[rowsin,colsin]
				} 
				if(r.overlap&c.overlap){
					X <- X - d0*u0%*%t(v0)
				}
				break
			}
			
		}
		if(stop){
			number <- k-1
			cat("no further stable bicluster, increase pcer!","\n")
			break
		}
		if(e==iter){
			number <- k-1
			cat("bicluster not coverged, increase number of iterations!","\n")
			break
		}
	}
	if(!number) number <- k
	niter[k] <- e
	Rspath[[k]] <- us[[2]]
	Cspath[[k]] <- vs[[2]]
	rows[rowsin] <- u0!=0
	cols[colsin] <- v0!=0
	Rows[[k]] <- rows
	Cols[[k]] <- cols
	ul[[k]] <- u0
	vl[[k]] <- v0
	dl[k] <- d0
	params <- list(steps = steps, pcer=pcer, iter=iter, ss.thr=ss.thr, size=size, weak=weak, r.overlap=r.overlap, c.overlap=c.overlap,
			r.negcorr=r.negcorr, c.negcorr=c.negcorr, nbiclust=nbiclust, merr=merr)
	RowxNumber=t(matrix(unlist(Rows),byrow=T,ncol=length(Rows[[1]])))
	NumberxCol=matrix(unlist(Cols),byrow=T,ncol=length(Cols[[1]]))
	Number <- number
	info <- list(ul=ul,vl=vl,dl=dl,Rspath=Rspath,Cspath=Cspath,niter=niter)
	return(BiclustResult(params,RowxNumber,NumberxCol,Number,info))
}

updateu <- function(X, v0, r.steps, r.err, ss.thr,size, weak, r.negcorr, mc.cores){
	p <- nrow(X)
	n <- ncol(X)
	ols <- X%*%v0
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	#lambdas <- seq(max(abs(ols)),0,length.out=n+1)
	lpath <- length(lambdas)    
	qumax <- sqrt((r.err)*((2*ss.thr)-1)*p)
	probpath <- matrix(nrow=p,ncol=lpath)
	if(r.negcorr){
		for(i in 1:lpath){
			ss.index <- sapply(1:r.steps,function(x){sample(1:n,n*size)})
			temp <- matrix(unlist((mclapply(1:r.steps,u.stepsnc,ss.index,v0,X,
												lambda=lambdas[i],weak,mc.cores=mc.cores))),nrow=p)
			qu <- mean(colSums(temp!=0))
			if(qu>qumax){ 
				break}
			probpath[,i] <- rowMeans(temp!=0)
		}
	}else{
		for(i in 1:lpath){
			ss.index <- sapply(1:r.steps,function(x){sample(1:n,n*size)})
			temp <- matrix(unlist((mclapply(1:r.steps,u.steps,ss.index,v0,X,
												lambda=lambdas[i],weak,mc.cores=mc.cores))),nrow=p)
			qu <- mean(colSums(temp!=0))
			if(qu>qumax){
				break}	
			probpath[,i] <- rowMeans(temp!=0)
		}
	}
	stable <- which(apply(probpath>=ss.thr,1,any))
	stop <- FALSE
	if(length(stable)==0) stop <- TRUE
	uc <- rep(0,p)
	uc[stable] <- (sign(ols)*(abs(ols)>=lambdas[i-1])*(abs(ols)-lambdas[i-1]))[stable]
	if(!r.negcorr) uc[which(sign(uc) != sign(sum(sign(uc)))) ] <- 0 
	return(list(uc,probpath,stop))
}

u.stepsnc <- function(ss,ss.index,v0,X,lambda,weak){
	ssX <- X[,ss.index[,ss]]
	ssv0 <- v0[ss.index[,ss]]
	ols <- ssX%*%ssv0
	delta <- lambda * runif(length(ols),weak,1) 
	uc <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	return(uc)
}	

u.steps <- function(ss,ss.index,v0,X,lambda,weak){
	ssX <- X[,ss.index[,ss]]
	ssv0 <- v0[ss.index[,ss]]
	ols <- ssX%*%ssv0
	delta <- lambda * runif(length(ols),weak,1) 
	uc <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	uc[which(sign(uc) != sign(sum(sign(uc)))) ] <- 0 
	return(uc)
}	

updatev <- function(X, u0, c.steps, c.err, ss.thr,size, weak, c.negcorr, mc.cores){
	p <- nrow(X)
	n <- ncol(X)
	ols <- t(X)%*%u0
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	#lambdas <- seq(max(abs(ols)),0,length.out=p+1)
	lpath <- length(lambdas)    
	qvmax <- sqrt((c.err)*((2*ss.thr)-1)*n)
	probpath <- matrix(nrow=n,ncol=lpath)
	if(c.negcorr){
		for(i in 1:lpath){
			ss.index <- sapply(1:c.steps,function(x){sample(1:p,p*size)})
			temp <- matrix(unlist((mclapply(1:c.steps,v.stepsnc,ss.index,u0,X,
												lambda=lambdas[i],weak,mc.cores=mc.cores))),nrow=n)
			qv <- mean(colSums(temp!=0))
			if(qv>qvmax){
				break
			}
			probpath[,i] <- rowMeans(temp!=0)
		}
	}else{
		for(i in 1:lpath){
			ss.index <- sapply(1:c.steps,function(x){sample(1:p,p*size)})
			temp <- matrix(unlist((mclapply(1:c.steps,v.steps,ss.index,u0,X,
												lambda=lambdas[i],weak,mc.cores=mc.cores))),nrow=n)
			qv <- mean(colSums(temp!=0))
			if(qv>qvmax){
				break
			}
			probpath[,i] <- rowMeans(temp!=0)
			}
	}
	stable <- which(apply(probpath>=ss.thr,1,any)) 
	stop <- FALSE
	if(length(stable)==0) stop <- TRUE
	vc <- rep(0,n)
	vc[stable] <- (sign(ols)*(abs(ols)>=lambdas[i-1])*(abs(ols)-lambdas[i-1]))[stable]
	if(!c.negcorr) vc[which(sign(vc) != sign(sum(sign(vc)))) ] <- 0 
	return(list(vc,probpath,stop))
}


v.stepsnc <- function(ss,ss.index,u0,X,lambda,weak){
	ssX <- X[ss.index[,ss],]
	ssu0 <- u0[ss.index[,ss]]
	ols <- t(ssX)%*%ssu0
	delta <- lambda * runif(length(ols),weak,1) 
	vc <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	return(vc)
}	

v.steps <- function(ss,ss.index,u0,X,lambda,weak){
	ssX <- X[ss.index[,ss],]
	ssu0 <- u0[ss.index[,ss]]
	ols <- t(ssX)%*%ssu0
	delta <- lambda * runif(length(ols),weak,1) 
	vc <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	vc[which(sign(vc) != sign(sum(sign(vc)))) ] <- 0 
	return(vc)
}	


