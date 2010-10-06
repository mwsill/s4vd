stabpath <- function(res,k){
	par(mfrow=c(1,2),omi=c(.1,.1,.3,.1)) 
	lRpath <- min(which(apply(res@info$Rspath[[k]],1,is.na)),ncol(res@info$Rspath[[k]]))
	cols <- ifelse(res@info$ul[[k]]==0,"black","red")
	if(k>res@Number) cols <- rep("black",length(res@info$ul))
	matplot(t(res@info$Rspath[[k]])[1:lRpath,],type="l",xlab=expression(paste(lambda[u])),
			ylab=paste("selection probability"),
			axes=F,lty=1,col=cols,ylim=c(0,1))
	axis(1,labels=FALSE)
	axis(2)
	title(paste(sum(res@info$ul[[k]]!=0)," rows"))
	abline(h=res@Parameters$ss.thr,col="red",lwd=3)
	text(x=lRpath/10,y=res@Parameters$ss.thr+0.05,labels=expression(paste(pi[thr])),col="red",cex=1.5)
	lCpath <- min(which(apply(res@info$Cspath[[k]],1,is.na)),ncol(res@info$Cspath[[k]]))
	cols <- ifelse(res@info$vl[[k]]==0,"black","red")
	if(k>res@Number) cols <- rep("black",length(res@info$vl))
	matplot(t(res@info$Cspath[[k]])[1:lCpath,],type="l",xlab=expression(paste(lambda[v])),
			ylab=paste("selection probability"),
			axes=F,lty=1,col=cols,ylim=c(0,1))
	axis(1,labels=FALSE)
	axis(2)
	title(paste(sum(res@info$vl[[k]]!=0)," columns"))
	abline(h=res@Parameters$ss.thr,col="red",lwd=3)
	text(x=lCpath/10,y=res@Parameters$ss.thr+0.05,labels=expression(paste(pi[thr])),col="red",cex=1.5)
	if(k>res@Number){
	title(paste("stability path bicluster",k," (not stable!)"),outer=T)	
	}else{
	title(paste("stability path bicluster",k),outer=T)
	}
	par(mfrow=c(1,1))
}

