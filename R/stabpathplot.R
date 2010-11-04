stabpath <- function(res,number){
	vc <- res@info[[number]][[1]]
	uc <- res@info[[number]][[2]]
	par(mfrow=c(1,2),omi=c(0.25, 0.25, 0.5, 0.25)) #c(bottom, left, top, right)
	n <- length(vc$qs)-1
	p <- length(uc$qs)-1
	pcer <- res@Parameters$pcer
	lv <- res@info[[number]][[1]]$l 
	lu <- res@info[[number]][[2]]$l
	thrv <- ((vc[[4]]^2/((n*pcer)*n))+1)/2
	thrv1 <- ((vc[[4]]^2/((n*pcer*2)*n))+1)/2
	thrv2 <- ((vc[[4]]^2/((n*0.5)*n))+1)/2
	redv <- which(vc[[2]][,lv]>thrv[lv])
	orangev <- which(vc[[2]][,lv]>thrv1[lv])
	greenv <- which(vc[[2]][,lv]>thrv2[lv])
	thru <- ((uc[[4]]^2/((p*pcer)*p))+1)/2
	thru1 <- ((uc[[4]]^2/((p*pcer*2)*p))+1)/2
	thru2 <- ((uc[[4]]^2/((p*0.5)*p))+1)/2
	redu <- which(uc[[2]][,lu]>thru[lu])
	orangeu <- which(uc[[2]][,lu]>thru1[lu])
	greenu <- which(uc[[2]][,lu]>thru2[lu])
	colsv <- rep("black",n)
	colsu <- rep("black",p)
	colsv[redv] <- "red"
	colsv[orangev] <- ifelse(colsv[orangev]=="black","orange",colsv[orangev])
	colsv[greenv] <- ifelse(colsv[greenv]=="black","green",colsv[greenv])
	colsu[redu] <- "red"
	colsu[orangeu] <- ifelse(colsu[orangeu]=="black","orange",colsu[orangeu])
	colsu[greenu] <- ifelse(colsu[greenu]=="black","green",colsu[greenu])
	matplot(t(vc[[2]]),type="l",col=colsv,lty=1,ylab="selection probability",xlab=expression(paste(lambda[v])),main="stability path columns",ylim=c(0,1))
	abline(v=lv,col="darkred",lwd=2)
	lines(thrv,col="darkred",lwd=3,lty=2)
	lines(thrv1,col="darkorange",lwd=3,lty=2)
	lines(thrv2,col="darkgreen",lwd=3,lty=2)
	legend(-(n*0.15), 1, c(paste("PCER ",pcer), paste("PCER ",pcer*2), "PCER 0.5"),text.col = c("darkred","darkorange","darkgreen"),bty="n")
	matplot(t(uc[[2]]),type="l",col=colsu,lty=1,ylab="selection probability",xlab=expression(paste(lambda[u])),main="stability path rows",ylim=c(0,1))
	abline(v=lu,col="darkred",lwd=2)
	lines(thru,col="darkred",lwd=3,lty=2)
	lines(thru1,col="darkorange",lwd=3,lty=2)
	lines(thru2,col="darkgreen",lwd=3,lty=2)
	legend(-(p*0.15), 1, c(paste("PCER ",pcer), paste("PCER ",pcer*2), "PCER 0.5"),text.col = c("darkred","darkorange","darkgreen"),bty="n")
	title(paste("Stability Paths Bicluster: ",number),outer=T)
}