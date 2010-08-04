s4vd <- function(
		X        		  # the p times n data matrix (p>>n)
		,r.steps=500      # number of subsamples used to calculate the relative selection frequencies for left singular vector coefficients that correspond to the rows of X
		,c.steps=500      # number of subsamples used to calculate the relative selection frequencies for right singular vector coefficients that correspond to the columns of X	
		,r.err=20         # expected number of falsly selected left singular vector coefficients 	
		,c.err=10         # expected number of falsly selected right singular vector coefficients 
		,iter=50	      # maximal number of iterations to fit a single bicluster	
		,ss.thr=0.6       # cutoff threshold (relative selection frequency) for the stability selection
		,size=0.5         # size of the subsamples for the stability selection
		,weak=0.2         # weakness parameter for the randomised lasso
		,r.overlap=TRUE   # allow bicluster to be column overlapping
		,c.overlap=TRUE   # allow bicluster to be row overlapping		
		,r.negcorr=TRUE   # allow for negative correlation of rows (genes) over columns (samples)  
		,c.negcorr=FALSE  # allow for negative correlation of columns (samples) over rows (genes)
		,nbiclust=50      # maximal number of biclusters
		,merr=0.01        # convergence parameter 
		,mc.cores=1       # number of cores for parallelization, computes the complete path
){
	startX <- X
	p <- nrow(X)
	n  <- ncol(X)
	rowsin <- rep(TRUE,nrow(X))
	colsin <- rep(TRUE,ncol(X))
	stop <- FALSE
	Rspath <- Cspath <- Rows <- Cols <- ul <- vl <- dl <- list()
	for(k in 1:nbiclust){
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
		cat("rows: ",sum(u0!=0)," columns: ",sum(v0!=0),"\n")
		if((length(u0)*size)<2|(length(v0)*size)<=2){
			cat("submatrix to small for resampling","\n")
			stop <- TRUE
			break
		}
		for(e in 1:iter){
			us <- updateu(X, v0, r.steps, r.err, ss.thr,size, weak, r.negcorr, mc.cores)
			u1 <- us[[1]]/sqrt(sum(us[[1]]^2))
			vs <- updatev(X, u1, c.steps, c.err, ss.thr,size, weak, c.negcorr, mc.cores)
			v1 <- vs[[1]]/sqrt(sum(vs[[1]]^2)) 
			ud <- sqrt(sum((u0-u1)^2))
			vd <- sqrt(sum((v0-v1)^2))
			r.in <- which(u1!=0)
			c.in <- which(v1!=0)
			d1 <- as.numeric(t(u0)%*%X%*%v0)
			cat("iter: ",e,"\n")
			cat("ud: ",ud," vd: ",vd,"\n")
			cat("rows: ",sum(u1!=0),"columns: ",sum(v1!=0),"\n")
			v0 <- v1
			u0 <- u1
			d0 <- d1
			if(vd < merr |ud < merr){
				cat("layer ",k," converged","\n")
				rows[rowsin] <- u0!=0
				cols[colsin] <- v0!=0
				Rspath[[k]] <- us[[2]]
				Cspath[[k]] <- vs[[2]]	
				Rows[[k]] <- rows
				Cols[[k]] <- cols
				ul[[k]] <- u0
				vl[[k]] <- v0
				dl[[k]] <- d0
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
		if(stop|e==iter)break
	}
	params <- list(r.steps = r.steps, c.steps=c.steps , r.err=r.err, c.err=c.err, iter=iter, ss.thr=ss.thr, size=size, weak=weak, r.overlap=r.overlap, c.overlap=c.overlap, r.negcorr=r.negcorr, c.negcorr=c.negcorr, nbiclust=nbiclust, merr=merr)
	RowxNumber=t(matrix(unlist(Rows),byrow=T,ncol=length(Rows[[1]])))
	NumberxCol=matrix(unlist(Cols),byrow=T,ncol=length(Cols[[1]]))
	Number <- length(Rows)
	info <- list(ul=ul,vl=vl,dl=dl,Rspath=Rspath,Cspath=Cspath)
	return(BiclustResult(params,RowxNumber,NumberxCol,Number,info))
}

updateu <- function(X, v0, r.steps, r.err, ss.thr,size, weak, r.negcorr, mc.cores){
	p <- nrow(X)
	n <- ncol(X)
	ols <- X%*%v0
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	lpath <- length(lambdas)    
	qumax <- sqrt((r.err)*((2*ss.thr)-1)*p)
	probpath <- matrix(nrow=p,ncol=lpath)
	if(r.negcorr){
		for(i in 1:lpath){
			ss.index <- sapply(1:r.steps,function(x){sample(1:n,n*size)})
			temp <- matrix(unlist((mclapply(1:r.steps,u.stepsnc,ss.index,v0,X,
												lambda=lambdas[i],weak,mc.cores=mc.cores))),nrow=p)
			qu <- mean(colSums(temp!=0))
			probpath[,i] <- rowMeans(temp!=0)
			if(qu>=qumax){ 
				break}	
		}
	}else{
		for(i in 1:lpath){
			ss.index <- sapply(1:r.steps,function(x){sample(1:n,n*size)})
			temp <- matrix(unlist((mclapply(1:r.steps,u.steps,ss.index,v0,X,
												lambda=lambdas[i],weak,mc.cores=mc.cores))),nrow=p)
			qu <- mean(colSums(temp!=0))
			probpath[,i] <- rowMeans(temp!=0)
			if(qu>=qumax){ 
				break}	
		}
	}
	stable <- which(apply(probpath>=ss.thr,1,any))
	uc <- sign(ols)*(abs(ols)>=lambdas[i-1])*(abs(ols)-lambdas[i-1])
	uc[-stable] <-0
	if(!r.negcorr) uc[which(sign(uc) != sign(sum(sign(uc)))) ] <- 0 
	return(list(uc,probpath))
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
	lpath <- length(lambdas)    
	qvmax <- sqrt((c.err)*((2*ss.thr)-1)*n)
	probpath <- matrix(nrow=n,ncol=lpath)
	if(c.negcorr){
		for(i in 1:lpath){
			ss.index <- sapply(1:c.steps,function(x){sample(1:p,p*size)})
			temp <- matrix(unlist((mclapply(1:c.steps,v.stepsnc,ss.index,u0,X,
												lambda=lambdas[i],weak,mc.cores=mc.cores))),nrow=n)
			qv <- mean(colSums(temp!=0))
			probpath[,i] <- rowMeans(temp!=0)
			if(qv>=qvmax){
				break
			}
		}
	}else{
		for(i in 1:lpath){
			ss.index <- sapply(1:c.steps,function(x){sample(1:p,p*size)})
			temp <- matrix(unlist((mclapply(1:c.steps,v.steps,ss.index,u0,X,
												lambda=lambdas[i],weak,mc.cores=mc.cores))),nrow=n)
			qv <- mean(colSums(temp!=0))
			probpath[,i] <- rowMeans(temp!=0)
			if(qv>=qvmax){
				break
			}
		}
	}
	stable <- which(apply(probpath>=ss.thr,1,any)) 
	if(length(stable)==0) stable <- which.max(apply(probpath,1,max)) 
	vc <- sign(ols)*(abs(ols)>=lambdas[i-1])*(abs(ols)-lambdas[i-1])
	vc[-stable] <-0
	if(!c.negcorr) vc[which(sign(vc) != sign(sum(sign(vc)))) ] <- 0 
	return(list(vc,probpath))
}


v.stepsnc <- function(ss,ss.index,u0,X,lambda,weak){
	ssX <- X[ss.index[,ss],]
	ssu0 <- u0[ss.index[,ss]]
	ols <- t(ssX)%*%ssu0
	delta <- lambda * runif(length(ols),1-weak,1) 
	vc <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	return(vc)
}	

v.steps <- function(ss,ss.index,u0,X,lambda,weak){
	ssX <- X[ss.index[,ss],]
	ssu0 <- u0[ss.index[,ss]]
	ols <- t(ssX)%*%ssu0
	delta <- lambda * runif(length(ols),1-weak,1) 
	vc <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	vc[which(sign(vc) != sign(sum(sign(vc)))) ] <- 0 
	return(vc)
}	