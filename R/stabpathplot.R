stabpath <- function(res,k){
	par(mfrow=c(1,2),omi=c(.1,.1,.3,.1)) 
	if(class(res@info$Rspath[[k]])=="numeric"){
		cols <- ifelse(res@info$ul[[k]]==0,"black","red")
		plot(res@info$Rspath[[k]],ylim=c(0,1),col=cols,yaxp=c(0.1,1,9),ylab=paste("selection probability"),xlab="rows")
		abline(h=0.5,col="black",lwd=1,lty=2)
		abline(h=res@info$thrs[[k]][2],col="red",lwd=3)
		title(paste(sum(res@info$ul[[k]]!=0)," rows"))
		text(x=nrow(res@RowxNumber)*0.9,y=res@info$thrs[[k]][2]+0.05,labels=expression(paste(pi[thr])),col="red",cex=1.5)
		cols <- ifelse(res@info$vl[[k]]==0,"black","red")
		plot(res@info$Cspath[[k]],ylim=c(0,1),col=cols,yaxp=c(0.1,1,9),ylab=paste("selection probability"),xlab="cols")
		abline(h=0.5,col="black",lwd=1,lty=2)
		abline(h=res@info$thrs[[k]][1],col="red",lwd=3)
		title(paste(sum(res@info$vl[[k]]!=0)," cols"))
		text(x=ncol(res@NumberxCol)*0.9,y=res@info$thrs[[k]][1]+0.05,labels=expression(paste(pi[thr])),col="red",cex=1.5)
		title(paste("pointwise selection probabilities for bicluster",k," PCER: ",res@Parameters$pcer),outer=T)
	}else{
	lRpath <- min(which(apply(res@info$Rspath[[k]],1,is.na)),ncol(res@info$Rspath[[k]]))
	cols <- ifelse(res@info$ul[[k]]==0,"black","red")
	if(k>res@Number) cols <- rep("black",length(res@info$ul))
	matplot(t(res@info$Rspath[[k]])[1:lRpath,],type="l",xlab=expression(paste(lambda[u])),
			ylab=paste("selection probability"),
			axes=T,lty=1,col=cols,ylim=c(0,1))
	#axis(1,labels=FALSE)
	#axis(2)
	title(paste(sum(res@info$ul[[k]]!=0)," rows"))
	abline(h=res@Parameters$ss.thr,col="red",lwd=3)
	text(x=lRpath/10,y=res@Parameters$ss.thr+0.05,labels=expression(paste(pi[thr])),col="red",cex=1.5)
	lCpath <- min(which(apply(res@info$Cspath[[k]],1,is.na)),ncol(res@info$Cspath[[k]]))
	cols <- ifelse(res@info$vl[[k]]==0,"black","red")
	if(k>res@Number) cols <- rep("black",length(res@info$vl))
	matplot(t(res@info$Cspath[[k]])[1:lCpath,],type="l",xlab=expression(paste(lambda[v])),
			ylab=paste("selection probability"),
			axes=T,lty=1,col=cols,ylim=c(0,1))
	#axis(1,labels=FALSE)
	#axis(2)
	title(paste(sum(res@info$vl[[k]]!=0)," columns"))
	abline(h=res@Parameters$ss.thr,col="red",lwd=3)
	text(x=lCpath/10,y=res@Parameters$ss.thr+0.05,labels=expression(paste(pi[thr])),col="red",cex=1.5)
	if(k>res@Number){
	title(paste("stability path bicluster",k," (not stable!)"),outer=T)	
	}else{
	title(paste("stability path bicluster",k," PCER: ",res@Parameters$pcer),outer=T)
	}
	}
	par(mfrow=c(1,1))
}

