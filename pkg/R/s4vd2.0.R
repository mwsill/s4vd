# TODO: Add comment
# 
# Author: martin
###############################################################################
s4vd <- function(
		X,
		steps = 100,
		pcer = 0.1,
		ss.thr = c(0.6,0.8),
		size = 0.5,
		gamm = 0,
		iter = 100,
		nbiclust = 20,
		merr = 0.001,
		cols.nn=FALSE,
		rows.nn=TRUE,
		row.overlap=TRUE,
		col.overlap=TRUE,
		fullpath=FALSE,
		pointwise=TRUE
){
	startX <- X
	p <- nrow(X)
	n <- ncol(X)
	rowsin <- rep(TRUE,p)	
	colsin <- rep(TRUE,n)
	stop <- FALSE
	info <- Rows <- Cols <- vc <- uc <- list()
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
		if(pointwise){
			for(i in 1:iter){
				vc <- updatev.pw(X,u0,pcer,ss.thr,steps,size,gamm,cols.nn,vc$l)
				v1 <- vc[[1]]/sqrt(sum(vc[[1]]^2)) 
				v1[is.na(v1)] <- 0
				uc <- updateu.pw(X,v1,pcer,ss.thr,steps,size,gamm,rows.nn,uc$l)
				u1 <- uc[[1]]/sqrt(sum(uc[[1]]^2)) 
				u1[is.na(u1)] <- 0
				if(vc[[3]]|uc[[3]]){
					stop <- TRUE
					break}
				ud <- sqrt(sum((u0-u1)^2))
				vd <- sqrt(sum((v0-v1)^2))
				cat("iter: ",i," rows: ",length(which(u1!=0))," cols: ",length(which(v1!=0))," merr: ",min(c(ud,vd)),"\n")
				r.in <- which(u1!=0)
				c.in <- which(v1!=0)
				v0 <- v1
				u0 <- u1
				if(min(c(vd,ud)) < merr)break
			}
		}else{
			for(i in 1:iter){
				vc <- updatev(X,u0,pcer,ss.thr,steps,size,gamm,cols.nn,fullpath)
				v1 <- vc[[1]]/sqrt(sum(vc[[1]]^2)) 
				v1[is.na(v1)] <- 0
				uc <- updateu(X,v1,pcer,ss.thr,steps,size,gamm,rows.nn,fullpath)
				u1 <- uc[[1]]/sqrt(sum(uc[[1]]^2)) 
				u1[is.na(u1)] <- 0
				if(vc[[3]]|uc[[3]]){
					stop <- TRUE
					break}
				ud <- sqrt(sum((u0-u1)^2))
				vd <- sqrt(sum((v0-v1)^2))
				cat("iter: ",i," rows: ",length(which(u1!=0))," cols: ",length(which(v1!=0))," merr: ",min(c(ud,vd)),"\n")
				r.in <- which(u1!=0)
				c.in <- which(v1!=0)
				v0 <- v1
				u0 <- u1
				if(min(c(vd,ud)) < merr)break
			}
		}
		rows[rowsin] <- u0!=0
		cols[colsin] <- v0!=0
		Rows[[k]] <- rows
		Cols[[k]] <- cols
		info[[k]] <- list(vc,uc)
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
	if(!stop) number <- k
	params <- list(steps = steps, pcer=pcer, iter=iter, ss.thr=ss.thr, size=size, gamm=gamm, row.overlap=row.overlap, col.overlap=col.overlap,
			rows.nn=rows.nn, cols.nn=cols.nn, nbiclust=nbiclust, merr=merr)
	RowxNumber=t(matrix(unlist(Rows),byrow=T,ncol=length(Rows[[1]])))
	NumberxCol=matrix(unlist(Cols),byrow=T,ncol=length(Cols[[1]]))
	Number <- number
	return(BiclustResult(params,RowxNumber,NumberxCol,Number,info))
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
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	set <- FALSE
	if(cols.nnc){
		for(l in 1:length(lambdas)){
			temp <- adaLasso.nn(t(X),u0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thrall[l] <- ((qs[l]^2/(err*n))+1)/2
			if(thrall[l]>=ss.thr[1]){
				if(!set){ ls <- l
				set <- TRUE}
				if(!fullpath)break
			}
		}
	}else{
		for(l in 1:length(lambdas)){
			temp <- adaLasso(t(X),u0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thrall[l] <- ((qs[l]^2/(err*n))+1)/2
			if(thrall[l]>=ss.thr[1]){
				if(!set){ ls <- l
				set <- TRUE}
				if(!fullpath)break
			}
		}
	}
	thr <- thrall[ls]
	if(thr > ss.thr[2]){
		while(pcer <= 0.5){
		pcer <- pcer + 0.01	
		thrall <- ((qs^2/((pcer*n)*n))+1)/2
		thr <- thrall[ls]
		if(thr < ss.thr[2])break
		}
	} 
	stable <- which(selprobpath[,ls]>=thr)
	if(length(stable)==0)stop <- TRUE
	vc <- rep(0,n)
	delta <- lambdas[ls]/(abs(ols)^gamm)  
	vc[stable] <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))[stable]
	return(list(vc=vc,selprobpath=selprobpath,stop=stop,qs=qs,thr=thr,l=ls))
}

#update u
updateu <- function(X,v0,pcer,ss.thr,steps,size,gamm,rows.nnc=FALSE,fullpath=FALSE){
	p <- nrow(X)
	err <- pcer*p
	ols <- X%*%v0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	set <- FALSE
	if(rows.nnc){
		for(l in 1:length(lambdas)){
			temp <- adaLasso.nn(X,v0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thrall[l] <- ((qs[l]^2/(err*p))+1)/2
			if(thrall[l]>=ss.thr[1]){
				if(!set){ ls <- l
						  set <- TRUE
					  	}
				if(!fullpath)break
			}
		}
	}else{
		for(l in 1:length(lambdas)){
			temp <- adaLasso(X,v0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thrall[l] <- ((qs[l]^2/(err*p))+1)/2
			if(thrall[l]>=ss.thr[1]){
				if(!set){ ls <- l
						  set <- TRUE
				}
				if(!fullpath)break
			}
		}
	}
	thr <- thrall[ls]
	if(thr > ss.thr[2]){
		while(pcer <= 0.5){
		pcer <- pcer + 0.01	
		thrall <- ((qs^2/((pcer*p)*p))+1)/2
		thr <- thrall[ls]
		if(thr < ss.thr[2])break
		}
	} 
	stable <- which(selprobpath[,ls]>=thr)
	if(length(stable)==0)stop <- TRUE
	uc <- rep(0,p)
	delta <- lambdas[ls]/(abs(ols)^gamm)  
	uc[stable] <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))[stable]
	return(list(uc=uc,selprobpath=selprobpath,stop=stop,qs=qs,thr=thr,l=ls))
}

#update v pointwise

updatev.pw <- function(X,u0,pcer,ss.thr,steps,size,gamm,cols.nnc=FALSE,l=NULL){
	n <- ncol(X)
	err <- pcer*n
	ols <- t(X)%*%u0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(is.null(l)) l <- which(lambdas==quantile(lambdas,0.5,type=1))[1]
	#search for a lambda
	if(cols.nnc){
		for(g in 1:(length(lambdas))){
			temp <- adaLasso(t(X),u0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thrall[l] <- ((qs[l]^2/(err*n))+1)/2
			#cat("qv:",qs[l],"piv: ",thrall[l], "lambda: ", l," PCERv:", err/n,"\n")
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){
				if(l == length(lambdas))break
				if(thrall[l+1]> ss.thr[2]){
					ls <- l+1
					break
				}
				l <- min(length(lambdas),l + ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0){ 
					l <- l-1
					if(l == 0)break
				} 
			}
			if(thrall[l] > ss.thr[2]){ 
				if(l == 0)break
				if(thrall[l-1]<ss.thr[1]&thrall[l-1]!=0){
					ls <- l
					break
				}
				l <- max(1,l - ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0){ 
					l <- l+1
					if(l == length(l))break
				} 
			}
		}
	}else{
		for(g in 1:(length(lambdas))){
			temp <- adaLasso.nn(t(X),u0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thrall[l] <- ((qs[l]^2/(err*n))+1)/2
			#cat("qv:",qs[l],"piv: ",thrall[l], "lambda: ", l," PCERv:", err/n,"\n")
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){
				if(l == length(lambdas))break
				if(thrall[l+1]> ss.thr[2]){
					ls <- l+1
					break
				}
				l <- min(length(lambdas),l + ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0){ 
					l <- l-1
					if(l == 0)break
				} 
			}
			if(thrall[l] > ss.thr[2]){ 
				if(l == 0)break
				if(thrall[l-1]<ss.thr[1]&thrall[l-1]!=0){
					ls <- l
					break
				}
				l <- max(1,l - ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0){ 
					l <- l+1
					if(l == length(l))break
				} 
			}
		}
	}
	thr <- thrall[ls]
	if(thr > ss.thr[2]){
		while(pcer <= 0.5){
			pcer <- pcer + 0.01	
			thrall <- ((qs^2/((pcer*n)*n))+1)/2
			thr <- thrall[ls]
			if(thr < ss.thr[2]){
				cat("PCERV: ",pcer,"\n")
				break}
		}
	} 
	stable <- which(selprobpath[,ls]>=thr)
	if(length(stable)==0)stop <- TRUE
	vc <- rep(0,n)
	delta <- lambdas[l]/(abs(ols)^gamm)  
	vc[stable] <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))[stable]
	return(list(vc=vc,selprobpath=selprobpath,stop=stop,qs=qs,thr=thr,l=ls,delta=delta))
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
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(is.null(l)) l <- which(lambdas==quantile(lambdas,0.5,type=1))[1]
	#search for a lambda
	if(rows.nnc){
		for(g in 1:(length(lambdas))){
			temp <- adaLasso(X,v0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thrall[l] <- ((qs[l]^2/(err*p))+1)/2
			#cat("qu:",qs[l],"piu: ",thrall[l], "lambda: ", l," PCERu:", err/p,"\n")
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){
				if(l == length(lambdas))break
				if(thrall[l+1]> ss.thr[2]){
					ls <- l+1
					break
				}
				l <- min(length(lambdas),l + ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0){ 
					l <- l-1
					if(l == 0)break
				} 
			}
			if(thrall[l] > ss.thr[2]){ 
				if(l == 0)break
				if(thrall[l-1]<ss.thr[1]&thrall[l-1]!=0){
					ls <- l
					break
				}
				l <- max(1,l - ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0){ 
					l <- l+1
					if(l == length(l))break
				} 
			}
		}
	}else{
		for(g in 1:(length(lambdas))){
			temp <- adaLasso.nn(X,v0,lambdas[l],steps,size,gamm)
			qs[l] <- mean(colSums(temp!=0))
			selprobpath[,l] <- rowMeans(temp!=0)
			thrall[l] <- ((qs[l]^2/(err*p))+1)/2
			#cat("qu:",qs[l],"piu: ",thrall[l], "lambda: ", l," PCERu:", err/p,"\n")
		if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
			ls <- l
			break
		} 
		if(thrall[l] < ss.thr[1]){
			if(l == length(lambdas))break
			if(thrall[l+1]> ss.thr[2]){
				ls <- l+1
				break
			}
			l <- min(length(lambdas),l + ceiling(length(lambdas)/(g+1)))
			while(thrall[l]!=0){ 
				l <- l-1
				if(l == 0)break
			} 
		}
		if(thrall[l] > ss.thr[2]){ 
			if(l == 0)break
			if(thrall[l-1]<ss.thr[1]&thrall[l-1]!=0){
				ls <- l
				break
			}
			l <- max(1,l - ceiling(length(lambdas)/(g+1)))
			while(thrall[l]!=0){ 
				l <- l+1
				if(l == length(l))break
			} 
		}
		}
	}
	thr <- thrall[l]
	if(thr > ss.thr[2]){
		while(pcer <= 0.5){
			pcer <- pcer + 0.01	
			thrall <- ((qs^2/((pcer*p)*p))+1)/2
			thr <- thrall[ls]
			if(thr < ss.thr[2]){
				cat("PCERU: ",pcer,"\n")
				break}
		}
	} 
	stable <- which(selprobpath[,ls]>=thr)
	if(length(stable)==0)stop <- TRUE
	uc <- rep(0,p)
	delta <- lambdas[l]/(abs(ols)^gamm)  
	uc[stable] <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))[stable]
	return(list(uc=uc,selprobpath=selprobpath,stop=stop,qs=qs,thr=thr,l=ls,delta=delta))
}

#adaptive Lasso 
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


