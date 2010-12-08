s4vd <- function(
		X,
		steps = 100,
		pcerv = 0.05,
		pceru = 0.05,
		ss.thr = c(0.6,0.65),
		size = 0.632,
		gamm = 0,
		iter = 100,
		nbiclust = 10,
		merr = 0.0001,
		cols.nc=FALSE,
		rows.nc=TRUE,
		row.overlap=TRUE,
		col.overlap=TRUE,
		row.min=4,
		col.min=4,
		pointwise=TRUE,
		start.iter=0,
		savepath=FALSE
){
	MYCALL<-match.call()
	startX <- X
	p.ini <- nrow(X)
	n.ini <- ncol(X)
	rowsin <- rep(TRUE,p.ini)	
	colsin <- rep(TRUE,n.ini)
	stop <- FALSE
	start <- TRUE
	info <- Rows <- Cols <- vc <- uc <- list()
	for(k in 1:nbiclust){
		gc()
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
		vc <- uc <- list()
		if((length(u0)*size)<2|(length(v0)*size)<=2){
			cat("submatrix to small for resampling","\n")
			number <- k-1
			stop <- TRUE
			break
		}
		if(pointwise){
			for(i in 1:iter){
				if(i > start.iter) start <- FALSE
				uc <- updateu.pw(X,v0,pceru,p.ini,ss.thr,steps,size,gamm,rows.nc,uc$l)
				u1 <- uc[[1]]/sqrt(sum(uc[[1]]^2)) 
				u1[is.na(u1)] <- 0
				vc <- updatev.pw(X,u1,pcerv,n.ini,ss.thr,steps,size,gamm,cols.nc,vc$l)
				v1 <- vc[[1]]/sqrt(sum(vc[[1]]^2)) 
				v1[is.na(v1)] <- 0
				if(uc[[3]] & i > start.iter){
					cat("rows not stable")
					stop <- TRUE
					break}
				if(vc[[3]] & i > start.iter){
					cat("columns not stable")
					stop <- TRUE
					break}
				ud <- sqrt(sum((u0-u1)^2))
				vd <- sqrt(sum((v0-v1)^2))
				cat("iter: ",i," rows: ",sum(uc$sp>=uc$thr)," cols: ",sum(vc$sp>=uc$thr)
						," merr: ",min(c(ud,vd)),"row.thr:",uc[[5]],"col.thr",vc[[5]],"\n")
				v0 <- v1
				u0 <- u1
				if(min(c(vd,ud)) < merr & i > start.iter)break
			}
		}else{
			for(i in 1:iter){
				if(i > start.iter) start <- FALSE
				uc <- updateu(X,v0,pceru,p.ini,ss.thr,steps,size,gamm,rows.nc,savepath)
				u1 <- uc[[1]]/sqrt(sum(uc[[1]]^2)) 
				u1[is.na(u1)] <- 0
				vc <- updatev(X,u1,pcerv,n.ini,ss.thr,steps,size,gamm,cols.nc,savepath)
				v1 <- vc[[1]]/sqrt(sum(vc[[1]]^2)) 
				v1[is.na(v1)] <- 0
				if(uc[[3]] & i > start.iter){
					cat("rows not stable")
					stop <- TRUE
					break}
				if(vc[[3]] & i > start.iter){
					cat("columns not stable")
					stop <- TRUE
					break}
				ud <- sqrt(sum((u0-u1)^2))
				vd <- sqrt(sum((v0-v1)^2))
				cat("iter: ",i," rows: ",sum(uc$sp>=uc$thr)," cols: ",sum(vc$sp>=uc$thr)
						," merr: ",min(c(ud,vd)),"row.thr:",uc[[5]],"col.thr",vc[[5]],"\n")
				v0 <- v1
				u0 <- u1
				if(min(c(vd,ud)) < merr & i > start.iter)break
			}
		}
		stableu <- uc$sp >= uc$thr
		stablev <- vc$sp >= vc$thr
		u0[!stableu] <- 0
		v0[!stablev] <- 0
		rows[rowsin] <- u0!=0
		cols[colsin] <- v0!=0
		Rows[[k]] <- rows
		Cols[[k]] <- cols
		if(sum(u0!=0)<row.min|sum(v0!=0)<col.min) stop <- TRUE
		if(!row.overlap){
			d0 <- as.numeric(t(u0)%*%X%*%v0)
			frobBC <- sqrt(sum((X - (d0*u0%*%t(v0)))[rowsin[rows],cols]^2))
			frobSVD <- sqrt(sum((X-SVD$d*SVD$u%*%t(SVD$v))[rowsin[rows],cols]^2))
			rowsin[rows] <- FALSE
			X <- startX[rowsin,colsin]
			info[[k]] <- list(vc,uc,frobenius=c(frobBC=frobBC,frobSVD=frobSVD))
		} 
		if(!col.overlap){
			d0 <- as.numeric(t(u0)%*%X%*%v0)
			frobBC <- sqrt(sum((X - (d0*u0%*%t(v0)))[rows,colsin[cols]]^2))
			frobSVD <- sqrt(sum((X-SVD$d*SVD$u%*%t(SVD$v))[rows,colsin[cols]]^2))
			colsin[cols] <- FALSE
			X <- startX[rowsin,colsin]
			info[[k]] <- list(vc,uc,frobenius=c(frobBC=frobBC,frobSVD=frobSVD))
		} 
		if(row.overlap&col.overlap){
			d0 <- as.numeric(t(u0)%*%X%*%v0)
			frobSVD <- sqrt(sum((X-(SVD$d*SVD$u%*%t(SVD$v)))[rows,cols]^2)) 
			X <- X - (d0*u0%*%t(v0))
			frobBC <- sqrt(sum(X[rows,cols]^2))
			info[[k]] <- list(vc,uc,frobenius=c(frobBC=frobBC,frobSVD=frobSVD))
		}
		if(stop){
			number <- k-1
			break
		}
		if(i==iter){
			number <- k-1
			cat("Fail to converge! Increase the number of iterations !","\n")
			gc()
			break
		}
	}
	if(!stop) number <- k
	params <- list(steps = steps, pcerv=pcerv, pceru=pceru, iter=iter, ss.thr=ss.thr, size=size, gamm=gamm, row.overlap=row.overlap, col.overlap=col.overlap,
			rows.nc=rows.nc, cols.nc=cols.nc, nbiclust=nbiclust, merr=merr, row.min=row.min, col.min=col.min, pointwise=pointwise, start.iter=start.iter, savepath=savepath, Call=MYCALL)  
	RowxNumber=t(matrix(unlist(Rows),byrow=T,ncol=length(Rows[[1]])))
	NumberxCol=matrix(unlist(Cols),byrow=T,ncol=length(Cols[[1]]))
	Number <- number
	return(BiclustResult(params,RowxNumber,NumberxCol,Number,info))
}


#update v
updatev <- function(X,u0,pcer,n.ini,ss.thr,steps,size,gamm,cols.nc=FALSE,savepath=FALSE){
	n <- ncol(X)
	err <- pcer*n.ini
	err*n
	ols <- t(X)%*%u0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	if(savepath) selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(cols.nc){
		for(l in 1:length(lambdas)){
			temp <- adaLasso.nc(t(X),u0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			if(savepath)selprobpath[,l] <- sp
			thrall[l] <- ((qs[l]^2/(err*n.ini))+1)/2
			if(thrall[l]>=ss.thr[1]){
				ls <- l
				break
			}
		}
	}else{
		for(l in 1:length(lambdas)){
			temp <- adaLasso(t(X),u0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			if(savepath)selprobpath[,l] <- sp
			thrall[l] <- ((qs[l]^2/(err*n.ini))+1)/2
			if(thrall[l]>=ss.thr[1]){
				ls <- l
				break
			}
		}
	}
	thr <- thrall[ls]
	if(thr > ss.thr[2]){
		while(pcer <= 0.5){
			pcer <- pcer + 0.01	
			thrall <- ((qs^2/((pcer*n.ini)*n.ini))+1)/2
			thr <- thrall[ls]
			if(thr < ss.thr[2])break
		}
	} 
	stable <- which(sp>=thr)
	if(length(stable)==0)stop <- TRUE
	vc <- rep(0,n)
	delta <- lambdas[ls]/(abs(ols)^gamm)  
	#if(savepath) sp <- selprobpath
	vc <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))
	if(savepath){
		return(list(vc=vc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,selprobpath=selprobpath))
	}else{
		return(list(vc=vc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls))
	}
}

#update u
updateu <- function(X,v0,pcer,p.ini,ss.thr,steps,size,gamm,rows.nc=FALSE,savepath=FALSE,start=FALSE){
	p <- nrow(X)
	err <- pcer*p
	ols <- X%*%v0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	if(savepath) selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(rows.nc){
		for(l in 1:length(lambdas)){
			temp <- adaLasso.nc(X,v0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			if(savepath)selprobpath[,l] <- sp
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			if(thrall[l]>=ss.thr[1]){
				ls <- l
				break
			}
		}
	}else{
		for(l in 1:length(lambdas)){
			temp <- adaLasso(X,v0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			if(savepath)selprobpath[,l] <- sp
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			if(thrall[l]>=ss.thr[1]){
				ls <- l
				break
			}
		}
	}
	thr <- thrall[ls]
	if(thr > ss.thr[2]){
		while(pcer <= 0.5){
			pcer <- pcer + 0.01	
			thrall <- ((qs^2/((pcer*p.ini)*p.ini))+1)/2
			thr <- thrall[ls]
			if(thr < ss.thr[2])break
		}
	} 
	stable <- which(sp>=thr)
	if(length(stable)==0)stop <- TRUE
	uc <- rep(0,p)
	delta <- lambdas[ls]/(abs(ols)^gamm)  
	uc <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))
	if(savepath){
		return(list(uc=uc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,selprobpath=selprobpath))
	}
	else{
		return(list(uc=uc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls))
	}	
}

#update v pointwise

updatev.pw <- function(X,u0,pcer,n.ini,ss.thr,steps,size,gamm,cols.nc=FALSE,l=NULL,start=FALSE){
	n <- ncol(X)
	err <- pcer*n
	ols <- t(X)%*%u0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(is.null(l)) l <- which(lambdas==quantile(lambdas,0.5,type=1))[1]
	#search for a lambda
	if(cols.nc){
		for(g in 1:(length(lambdas))){
			temp <- adaLasso.nc(t(X),u0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*n.ini))+1)/2
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){
				if(l == length(lambdas))break
				if(thrall[l+1]> ss.thr[2]){
					ls <- l+1
					temp <- adaLasso.nc(t(X),u0,lambdas[ls],steps,size,gamm)
					t <- temp!=0
					qs[ls] <- mean(colSums(t))
					sp <- rowMeans(t)
					thrall[ls] <- ((qs[ls]^2/(err*n.ini))+1)/2
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
			temp <- adaLasso(t(X),u0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*n.ini))+1)/2
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){
				if(l == length(lambdas))break
				if(thrall[l+1]> ss.thr[2]){
					ls <- l+1
					temp <- adaLasso(t(X),u0,lambdas[ls],steps,size,gamm)
					t <- temp!=0
					qs[ls] <- mean(colSums(t))
					sp <- rowMeans(t)
					thrall[ls] <- ((qs[l]^2/(err*n.ini))+1)/2
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
			thrall <- ((qs^2/((pcer*n.ini)*n.ini))+1)/2
			thr <- thrall[ls]
			if(thr < ss.thr[2])	break
		}
	}
	thr <- ((qs[ls]^2/((pcer*n.ini)*n.ini))+1)/2
	stable <- which(sp>=thr)
	if(length(stable)==0)stop <- TRUE
	vc <- rep(0,n)
	delta <- lambdas[l]/(abs(ols)^gamm)  
	vc <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))
	return(list(vc=vc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,delta=delta))
}


#update u pointwise
updateu.pw <- function(X,v0,pcer,p.ini,ss.thr,steps,size,gamm,rows.nc=FALSE,l=NULL,start=FALSE){
	p <- nrow(X)
	#err <- pcer*p.ini
	err <- pcer*p
	ols <- X%*%v0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	qs <- numeric(length(lambdas)) 
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(is.null(l)) l <- which(lambdas==quantile(lambdas,0.5,type=1))[1]
	#search for a lambda
	if(rows.nc){
		for(g in 1:(length(lambdas))){
			temp <- adaLasso.nc(X,v0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){
				if(l == length(lambdas))break
				if(thrall[l+1]> ss.thr[2]){
					ls <- l+1
					temp <- adaLasso.nc(X,v0,lambdas[ls],steps,size,gamm)
					t <- temp!=0
					qs[ls] <- mean(colSums(t))
					sp <- rowMeans(t)
					thrall[ls] <- ((qs[ls]^2/(err*p.ini))+1)/2
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
			temp <- adaLasso(X,v0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){
				if(l == length(lambdas))break
				if(thrall[l+1]> ss.thr[2]){
					ls <- l+1
					temp <- adaLasso(X,v0,lambdas[ls],steps,size,gamm)
					t <- temp!=0
					qs[ls] <- mean(colSums(t))
					sp <- rowMeans(t)
					thrall[ls] <- ((qs[ls]^2/(err*p.ini))+1)/2
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
			thrall <- ((qs^2/((pcer*p.ini)*p.ini))+1)/2
			thr <- thrall[ls]
			if(thr < ss.thr[2])break
		}
	}
	thr <- ((qs[ls]^2/((pcer*p.ini)*p.ini))+1)/2
	stable <- which(sp>=thr)
	if(length(stable)==0)stop <- TRUE
	uc <- rep(0,p)
	delta <- lambdas[l]/(abs(ols)^gamm)  
	uc <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))
	return(list(uc=uc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,delta=delta))
}

#adaptive Lasso 
adaLasso.nc <- function(X,b,lambda,steps,size,gamm=0){
	subsets <- sapply(1:steps,function(x){sample(1:length(b),length(b)*size)})
	res <- sapply(1:steps,adaLassoSteps.nc,subsets,X,b,lambda,gamm)
	return(res)
}
adaLassoSteps.nc <- function(index,subsets,X,b,lambda,gamm){
	ols <- X[,subsets[,index]]%*%b[subsets[,index]]
	delta <- lambda/(abs(ols)^gamm)                        
	ols <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	ols[is.na(ols)] <- 0
	return(ols)
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
	mostof <- sign(sum(sign(ols)))
	if(mostof==0) mostof <- 1
	ols[which(sign(ols) != mostof) ] <- 0
	ols[is.na(ols)] <- 0
	return(ols)
}
