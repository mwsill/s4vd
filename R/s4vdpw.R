s4vdpw <- function(
		X        		    # data matrix 
		,steps=200          # number of subsamples used to calculate the relative selection frequencies for left singular vector coefficients that correspond to the rows of X
		,pcer=0.1           # desired per comparison error rate  	
		,pcer.max=2*pcer    #  
		,iter=100	        # maximal number of iterations to fit a single bicluster	
		,ss.thr=c(0.6,0.7)  # desired area of cutoff threshold (relative selection frequency) for the stability selection
		,size=0.5           # size of the subsamples for the stability selection
		,gamm=1             # gamma parmaeter for the adaptive lasso 0 == lasso  
		,r.overlap=TRUE     # allow bicluster to be column overlapping
		,c.overlap=TRUE     # allow bicluster to be row overlapping		
		,r.negcorr=TRUE     # allow for negative correlation of rows (genes) over columns (samples)  
		,c.negcorr=FALSE     # allow for negative correlation of columns (samples) over rows (genes)
		,nbiclust=20        # maximal number of biclusters
		,merr=0.05          # convergence parameter 
		,mc.cores=1         # number of cores for parallelization   
		,verbose=0          # verbose level 0 - 1 - 2
){
	startX <- X
	number <- FALSE
	us <- list(NULL,NULL,FALSE)
	p <- nrow(X)
	n  <- ncol(X)
	rowsin <- rep(TRUE,nrow(X))
	colsin <- rep(TRUE,ncol(X))
	stop <- FALSE
	Rspath <- Cspath <- Rows <- Cols <- ul <- vl <- thrs <- list()
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
		lv <- NULL
		lu <- NULL
		for(e in 1:iter){
			cat(".")
			pcerv <- pceru <- pcer
			#update v
			vs <- updatev.pw(X, u0, steps, pcer,pcer.max,ss.thr, size, gamm, c.negcorr, lv,G=G,verbose)
			lv <- vs[[4]]
			if(vs[[3]]|all(vs[[1]]==0)){
				stop <- TRUE
				if(pcer.max>pcer){
					while(TRUE){
						pcerv <- pcerv +0.01
						if(verbose>1)cat("PCERV: ", pcerv,"\n")
						vs <- updatev.pw(X, u0, steps, pcerv,pcer.max,ss.thr, size, gamm, c.negcorr,lv,G=G,verbose)
						lv <- vs[[4]]
						if(!vs[[3]]|pcerv>=pcer.max){
							if(!vs[[3]])stop <- FALSE
							break
						}
					}
				}
				if(stop)break
			}
			v1 <- vs[[1]]/sqrt(sum(vs[[1]]^2)) 
			v1[is.na(v1)] <- 0
			lv <- vs[[4]]
			
			# update u
			us <- updateu.pw(X, v1, steps, pcer,pcer.max,ss.thr,size, gamm, r.negcorr, lu,G=G,verbose)
			lu <- us[[4]]
			if(us[[3]]|all(us[[1]]==0)){
				stop <- TRUE
				if(pcer.max>pcer){
					while(TRUE){
						pceru <- pceru +0.01
						if(verbose>1)cat("PCERU: ", pceru,"\n")
						us <- updateu.pw(X, v1, steps, pceru,pcer.max,ss.thr,size, gamm, r.negcorr, lu,G=G,verbose)
						lu <- us[[4]]
						if(!us[[3]]|pceru>=pcer.max){
							if(!us[[3]])stop <- FALSE
							break
						}
					}
				}
				if(stop)break
			}
			u1 <- us[[1]]/sqrt(sum(us[[1]]^2))
			u1[is.na(u1)] <- 0
			
			
			ud <- sqrt(sum((u0-u1)^2))
			vd <- sqrt(sum((v0-v1)^2))
			if(verbose > 0) cat("iter: ",e," rows: ",length(which(u1!=0))," cols: ",length(which(v1!=0))," merr: ",min(c(ud,vd))," PCERu: ",pceru," PCERv: ",pcerv,"\n")
			r.in <- which(u1!=0)
			c.in <- which(v1!=0)
			d1 <- as.numeric(t(u0)%*%X%*%v0)
			v0 <- v1
			u0 <- u1
			d0 <- d1
			if(min(c(vd,ud)) < merr & e > 1 & pceru == pcer & pcerv == pcer){
				cat("> nrows:",sum(u0!=0),"ncols: ",sum(v0!=0),"\n")
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
				thrs[[k]] <- c(vs[[5]],us[[5]])
				if(!r.overlap){
					rowsin[rows] <- FALSE
					X <- startX[rowsin,colsin]
				} 
				if(!c.overlap){
					colsin[cols] <- FALSE
					X <- startX[rowsin,colsin]
				} 
				if(r.overlap&c.overlap){
					X <- round(X - d0*u0%*%t(v0),digits=10)
				}
				break
			}
			
		}
		if(stop){
			number <- k-1
			cat("no further stable bicluster, increase PCER!","\n")
			break
		}
		if(e==iter){
			number <- k-1
			cat("bicluster not coverged, increase merr or the number of iterations!","\n")
			break
		}
	}
	if(!number) number <- k
	niter[k] <- e
	Rspath[[k]] <- us[[2]]
	Cspath[[k]] <- vs[[2]]
	Rows[[k]] <- rows
	Cols[[k]] <- cols
	ul[[k]] <- u0
	vl[[k]] <- v0
	dl[k] <- d0
	params <- list(steps = steps, pcer=pcer, iter=iter, ss.thr=ss.thr, size=size, gamm=gamm, r.overlap=r.overlap, c.overlap=c.overlap,
			r.negcorr=r.negcorr, c.negcorr=c.negcorr, nbiclust=nbiclust, merr=merr)
	RowxNumber=t(matrix(unlist(Rows),byrow=T,ncol=length(Rows[[1]])))
	NumberxCol=matrix(unlist(Cols),byrow=T,ncol=length(Cols[[1]]))
	Number <- number
	info <- list(ul=ul,vl=vl,dl=dl,Rspath=Rspath,Cspath=Cspath,niter=niter,thrs=thrs)
	return(BiclustResult(params,RowxNumber,NumberxCol,Number,info))
}

updateu.pw <- function(X, v0, r.steps, pcer,pcer.max, ss.thr, size, gamm, r.negcorr, l=NULL, G=10,verbose){
	p <- nrow(X)
	n <- ncol(X)
	err <- pcer*p
	ols <- X%*%v0
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	if(is.null(l)) l <- which(lambdas==quantile(lambdas,0.5,type=1))[1]
	for(g in 1:(length(lambdas))){
		lambda <- lambdas[l]
		ss.index <- sapply(1:r.steps,function(x){sample(1:n,n*size)})
		if(r.negcorr){
			temp <- matrix(unlist((mclapply(1:r.steps,u.stepsnc,ss.index,v0,X,
												lambda=lambda,gamm))),nrow=p)
			qu <- mean(colSums(temp!=0))
			thr <- ((qu^2/(err*p))+1)/2
		}
		else{
			temp <- matrix(unlist((mclapply(1:r.steps,u.steps,ss.index,v0,X,
												lambda=lambda,gamm))),nrow=p)
			qu <- mean(colSums(temp!=0))
			thr <- ((qu^2/(err*p))+1)/2
		}
		if(verbose>1) cat("qu: ",qu,"piu: ",thr,"lambda: ", l,"\n")
		if(l == length(lambdas) & thr < ss.thr[1]){ 
			while(pcer > 0.02){
				pcer <- pcer - 0.01
				err <- pcer * n
				thr <- ((qu^2/(err*n))+1)/2
				if(thr >= ss.thr[1] & thr <= ss.thr[2])break 
			}
			if(verbose>1) cat("qu:",qu,"piu: ",thr, "lambda: ", l," PCERU:", pcer,"\n")
		}
		if(l==1 & thr > ss.thr[2]){ 
			while(pcer < pcer.max){
				pcer <- pcer + 0.01
				err <- pcer * n
				thr <- ((qu^2/(err*n))+1)/2
				if(thr >= ss.thr[1] & thr <= ss.thr[2])break 
			}
			if(verbose>1) cat("qu:",qu,"piu: ",thr, "lambda: ", l," PCERU:", pcer,"\n")
		}
		#if(thr==0.5 & l==length(lambdas))break
		#if(thr==0.5) l  <- length(lambdas)
		if(thr >= ss.thr[1] & thr <= ss.thr[2])break 
		if(thr < ss.thr[1]) l <- min(length(lambdas),l + ceiling(length(lambdas)/(g+1))) 
		if(thr > ss.thr[2]) l <- max(1,l - ceiling(length(lambdas)/(g+1))) 
	}
	stable <- which(rowMeans(temp!=0)>=thr)
	stop <- FALSE
	if(length(stable)==0)stop <- TRUE
	uc <- rep(0,p)
	delta <- lambda/(abs(ols)^gamm)  
	uc[stable] <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))[stable]
	if(!r.negcorr) uc[which(sign(uc) != names(which.max(table(sign(uc))))) ] <- 0
	probpath <- rowMeans(temp!=0)
	return(list(uc,probpath,stop,l,thr))
}



u.stepsnc <- function(ss,ss.index,v0,X,lambda,gamm){
	ssX <- X[,ss.index[,ss]]
	ssv0 <- v0[ss.index[,ss]]
	ols <- ssX%*%ssv0
	delta <- lambda/(abs(ols)^gamm)                      #adaptive lasso 
	#delta <- lambda * runif(length(ols),gamm,1)       #randomised lasso
	uc <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	uc[is.na(uc)] <- 0
	return(uc)
}	

u.steps <- function(ss,ss.index,v0,X,lambda,gamm){
	ssX <- X[,ss.index[,ss]]
	ssv0 <- v0[ss.index[,ss]]
	ols <- ssX%*%ssv0
	delta <- lambda/(abs(ols)^gamm)                      #adaptive lasso 
	#delta <- lambda * runif(length(ols),gamm,1)       #randomised lasso
	uc <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	uc[which(sign(uc) != names(which.max(table(sign(uc))))) ] <- 0
	uc[is.na(uc)] <- 0
	return(uc)
}	

updatev.pw <- function(X, u0, c.steps, pcer,pcer.max, ss.thr,size, gamm, c.negcorr, l=NULL, G=10, verbose){
	p <- nrow(X)
	n <- ncol(X)
	err <- pcer * n
	ols <- t(X)%*%u0
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	if(is.null(l)) l <- which(lambdas==quantile(lambdas,0.3,type=1))[1]
	for(g in 1:(length(lambdas))){
		lambda <- lambdas[l]
		ss.index <- sapply(1:c.steps,function(x){sample(1:p,p*size)})
		if(c.negcorr){
			temp <- matrix(unlist((mclapply(1:c.steps,v.stepsnc,ss.index,u0,X,
												lambda=lambda,gamm))),nrow=n)
			qv <- mean(colSums(temp!=0))
			thr <- ((qv^2/(err*n))+1)/2
		}
		else{
			temp <- matrix(unlist((mclapply(1:c.steps,v.steps,ss.index,u0,X,
												lambda=lambda,gamm))),nrow=n)
			qv <- mean(colSums(temp!=0))
			thr <- ((qv^2/(err*n))+1)/2
		}
		if(verbose>1) cat("qv:",qv,"piv: ",thr,"lambda: ", l,"\n")
		if(l == length(lambdas) & thr < ss.thr[1]){ 
			while(pcer > 0.02){
				pcer <- pcer - 0.01
				err <- pcer * p
				thr <- ((qv^2/(err*p))+1)/2
				if(thr >= ss.thr[1] & thr <= ss.thr[2])break 
			}
			if(verbose>1) cat("qv:",qv,"piv: ",thr,"lambda: ",l," PCERV:", pcer,"\n")
		}
		if( l==1 & thr > ss.thr[2]){ 
			while(pcer < pcer.max){
				pcer <- pcer + 0.01
				err <- pcer * p
				thr <- ((qv^2/(err*p))+1)/2
				if(thr >= ss.thr[1] & thr <= ss.thr[2])break 
			}
			if(verbose>1) cat("qv:",qv,"piv: ",thr,"lambda: ",l," PCERV:", pcer,"\n")
		}
		#	if(thr==0.5 & l==length(lambdas))break
	    #	if(thr==0.5) l  <- length(lambdas)
		if(thr >= ss.thr[1] & thr <= ss.thr[2])break 
		if(thr < ss.thr[1]) l <- min(length(lambdas),l + ceiling(length(lambdas)/(g+1))) 
		if(thr > ss.thr[2]) l <- max(1,l - ceiling(length(lambdas)/(g+1))) 
	}
	stable <- which(rowMeans(temp!=0)>=thr)
	stop <- FALSE
	if(length(stable)==0)stop <- TRUE
	vc <- rep(0,n)
	delta <- lambda/(abs(ols)^gamm) 
	vc[stable] <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))[stable]
	if(!c.negcorr) vc[which(sign(vc) != names(which.max(table(sign(vc))))) ] <- 0 
	probpath <- rowMeans(temp!=0)
	return(list(vc,probpath,stop,l,thr))
}

v.stepsnc <- function(ss,ss.index,u0,X,lambda,gamm){
	ssX <- X[ss.index[,ss],]
	ssu0 <- u0[ss.index[,ss]]
	ols <- t(ssX)%*%ssu0
	delta <- lambda/(abs(ols)^gamm)                      #adaptive lasso 
	#delta <- lambda * runif(length(ols),gamm,1)       #randomised lasso
	vc <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	vc[is.na(vc)] <- 0
	return(vc)
}	

v.steps <- function(ss,ss.index,u0,X,lambda,gamm){
	ssX <- X[ss.index[,ss],]
	ssu0 <- u0[ss.index[,ss]]
	ols <- t(ssX)%*%ssu0
	delta <- lambda/(abs(ols)^gamm)                      #adaptive lasso 
	#delta <- lambda * runif(length(ols),gamm,1)       #randomised lasso
	vc <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	vc[which(sign(vc) != names(which.max(table(sign(vc))))) ] <- 0 
	vc[is.na(vc)] <- 0
	return(vc)
}	




