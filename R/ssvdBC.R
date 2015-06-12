ssvdBC <- 
function(X, K=10, threu = 1, threv = 1,
         gamu = 0, gamv =0 , merr = 1e-4, niter = 100)
{
    MYCALL <- match.call()
    res <- list()
    RowxNumber <- matrix(nrow=nrow(X),ncol=K)
    NumberxCol <- matrix(ncol=ncol(X),nrow=K)
    for(k in 1:K){ 
## Two notes:
## 1. The original computation here is incorrect:
##      res[[k]]  = ssvd(X,threu = 1, threv = 1, gamu = 0, gamv =0,  u0 = svd(X)$u[,k], v0 = svd(X)$v[,k], merr = 10^(-4), niter = 100)
## Instead, it should have been:
##      res[[k]]  = ssvd(X,threu = 1, threv = 1, gamu = 0, gamv =0,  u0 = svd(X)$u[,1], v0 = svd(X)$v[,1], merr = 10^(-4), niter = 100)
## according to the deflation method defined in the paper of Lee, Shen, et. al.

# 2. We leave off the initial vectors to let our modified ssvd function compute
#    them using a more efficient partial SVD, and pass through the user-supplied
#    parameters that seem to be neglected in the original:
      res[[k]]  <- ssvd(X,threu=threu, threv=threv,
                                  gamu=gamu,gamv=gamv,merr=merr,niter=niter)
      if(res[[k]]$stop){
        K <- k-1
        break
      }
      RowxNumber[,k] <- res[[k]][[1]]!=0
      NumberxCol[k,] <- res[[k]][[2]]!=0
      d <- as.numeric( (t(res[[k]][[1]])%*%X%*%res[[k]][[2]]) [])
      res[[k]][[4]] <- d
#     X <- X - (d*res[[k]][[1]]%*%t(res[[k]][[2]]))
      X <- X - d*tcrossprod(res[[k]][[1]], res[[k]][[2]])
    }
    params <- list(K=K,threu=threu,threv=threv,gamu=gamv,merr=merr,niter=niter,Call=MYCALL)
    Number <- K
    info <- list(res=res)
    BiclustResult(params,RowxNumber,NumberxCol,Number,info)
}
