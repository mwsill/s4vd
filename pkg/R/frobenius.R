frobenius <- function(res){
	#if(!class(res@Parameters$Method)=="BCs4vd"){
	#	stop("object is not of class BCs4vd")
	#}
	frobBC <- numeric()
	frobSVD <- numeric()
	frobBCsvd <- numeric()
	for(k in 1:res@Number){ 
		 frobBC[k] <- res@info[[k]][[3]][1]
		 frobSVD[k] <- res@info[[k]][[3]][2]
		 frobBCsvd[k] <- res@info[[k]][[3]][3]
	 }
	 frob <- rbind(frobBC,frobSVD,frobBCsvd)
	 rownames(frob) <- c("Frobenius S4VD:","Frobenius SVD:","Frobenius SVD new:")
	 colnames(frob) <- paste("BC", 1:res@Number)
	 print(frob)
	 invisible(frob)
}