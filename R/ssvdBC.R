ssvdBC <- function(X,K=10,threu = 1, threv = 1, gamu = 0, gamv =0 , merr = 10^(-4), niter = 100){
		MYCALL <- match.call()
		res <- list()
		RowxNumber <- matrix(nrow=nrow(X),ncol=K)
		NumberxCol <- matrix(ncol=ncol(X),nrow=K)
		for(k in 1:K){ 
			res[[k]]	<- ssvd(X,threu = 1, threv = 1, gamu = 0, gamv =0,  u0 = svd(X)$u[,k], v0 = svd(X)$v[,k], merr = 10^(-4), niter = 100)
			if(res[[k]]$stop){
				K <- k-1
				break
			}
			RowxNumber[,k] <- res[[k]][[1]]!=0
			NumberxCol[k,] <- res[[k]][[2]]!=0
			d <- as.integer(t(res[[k]][[1]])%*%X%*%res[[k]][[2]])
			res[[k]][[4]] <- d
			X <- X - (d*res[[k]][[1]]%*%t(res[[k]][[2]]))
		}
		params <- list(K=K,threu=threu,threv=threv,gamu=gamv,merr=merr,niter=niter,Call=MYCALL)
		Number <- K
		info <- list(res=res)
		return(BiclustResult(params,RowxNumber,NumberxCol,Number,info))
}
