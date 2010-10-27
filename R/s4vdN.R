# TODO: Add comment
# 
# Author: martin
###############################################################################

dat <- datagen()
X <- dat[[1]]
steps <- 1000
size <- 0.5
pcer <- 0.1
gamm <- 0
SVD <- svd(X,nu=1,nv=1)
v0 <- SVD$v
u0 <- SVD$u
d0 <- SVD$d

p <- nrow(X)
n  <- ncol(X)

#update u
err <- pcer*p
ols <- X%*%v0
lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
#lambdas <- seq(max(abs(ols)),0,length.out = p+1)
selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
qus <- numeric(length(lambdas)) 
for(l in 1:length(lambdas)){
	cat(l,"\n")
	temp <- adaLasso(X,v0,lambdas[l],steps,size,gamm)
	qus[l] <- mean(colSums(temp!=0))
	selprobpath[,l] <- rowMeans(temp!=0)
}
thr3 <- ((qus^2/(10*p))+1)/2
thr2 <- ((qus^2/(20*p))+1)/2
thr1 <- ((qus^2/(30*p))+1)/2
thr <- ((qus^2/(50*p))+1)/2
matplot(t(selprobpath),type="l",col="black",lty=1)
lines(thr1,col="green",lwd=3)
lines(thr2,col="orange",lwd=3)
lines(thr3,col="red",lwd=3)
lines(thr,col="blue",lwd=3)

#updatev
ols <- t(X)%*%u0
lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
#lambdas <- seq(max(abs(ols)),0,length.out = p+1)
l <- 50

temp <- adaLasso(t(X),u0,lambda[l],steps,size,gamm)
qv <- mean(colSums(temp!=0))
thr <- ((qv^2/(err*n))+1)/2


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


t(X)%*%u0