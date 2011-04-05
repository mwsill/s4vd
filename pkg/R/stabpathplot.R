stabpath <- function(res,number){
	#if(!class(res@Parameters$Method)=="BCs4vd"){
	#	stop("object is not of class BCs4vd")
	#}
	if(res@Parameters$savepath&!res@Parameters$pointwise){
		vc <- res@info[[number]][[1]]
		uc <- res@info[[number]][[2]]
		par(mfrow=c(1,2),omi=c(0.25, 0.25, 0.5, 0.25)) #c(bottom, left, top, right)
		n <- length(vc$qs)-1
		p <- length(uc$qs)-1
		pcerv <- res@Parameters$pcerv
		lv <- res@info[[number]][[1]]$l 
		lu <- res@info[[number]][[2]]$l
		thrv <- ((vc[[4]]^2/((n*pcerv)*n))+1)/2
		redv <- which(vc[[7]][,lv]>thrv[lv])
		pceru <- res@Parameters$pceru
		thru <- ((uc[[4]]^2/((p*pceru)*p))+1)/2
		redu <- which(uc[[7]][,lu]>thru[lu])
		colsv <- rep("black",n)
		colsu <- rep("black",p)
		colsv[redv] <- "red"
		colsu[redu] <- "red"
		matplot(t(vc[[7]]),type="l",col=colsv,lty=1,ylab="selection probability",xlab=expression(paste(lambda[v])),main="stability path columns",ylim=c(0,1))
		abline(v=lv,col="darkred",lwd=2)
		abline(h=thrv[lv],col="darkred")
		lines(thrv,col="darkred",lwd=3,lty=2)
		legend(-(n*0.15), 1.05, c(paste("PCER ",pcerv)),text.col = c("darkred"),bty="n")
		matplot(t(uc[[7]]),type="l",col=colsu,lty=1,ylab="selection probability",xlab=expression(paste(lambda[u])),main="stability path rows",ylim=c(0,1))
		abline(v=lu,col="darkred",lwd=2)
		abline(h=thru[lu],col="darkred")
		lines(thru,col="darkred",lwd=3,lty=2)
		legend(-(p*0.15), 1.05, c(paste("PCER ",pceru)),text.col = c("darkred"),bty="n")
		title(paste("Stability Paths Bicluster: ",number),outer=T)
	}else{
		vc <- res@info[[number]][[1]]
		uc <- res@info[[number]][[2]]
		par(mfrow=c(1,2),omi=c(0.25, 0.25, 0.5, 0.25)) #c(bottom, left, top, right)
		n <- length(vc$qs)-1
		p <- length(uc$qs)-1
		pcerv <- res@Parameters$pcerv
		pceru <- res@Parameters$pceru
		lv <- res@info[[number]][[1]]$l 
		lu <- res@info[[number]][[2]]$l
		redv <- which(vc[[2]]>=res@info[[number]][[1]]$thr)
		redu <- which(uc[[2]]>=res@info[[number]][[2]]$thr)
		colsv <- rep("black",n)
		colsu <- rep("black",p)
		colsv[redv] <- "red"
		colsu[redu] <- "red"
		plot(vc[[2]],col=colsv,cex=.5,ylim=c(0,1),ylab="selection probability",main="stability path columns",xlab="columns")
		abline(h=res@info[[number]][[1]]$thr,col="darkred")
		legend(-(n*0.15), 1.05, c(paste("PCER ",pcerv)),text.col = c("darkred"),bty="n")
		plot(uc[[2]],col=colsu,cex=.5,ylim=c(0,1),ylab="selection probability",main="stability path rows",xlab="rows")
		abline(h=res@info[[number]][[2]]$thr,col="darkred")
		legend(-(p*0.15), 1.05, c(paste("PCER ",pceru)),text.col = c("darkred"),bty="n")
	}
}