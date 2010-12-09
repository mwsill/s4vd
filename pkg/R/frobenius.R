frobenius <- function(res){
	#if(!class(res@Parameters$Method)=="BCs4vd"){
	#	stop("object is not of class BCs4vd")
	#}
	frobBC <- numeric()
	frobSVD <- numeric()		
	for(k in 1:res@Number){ 
		 frobBC[k] <- res@info[[k]][[3]][1]
		 frobSVD[k] <- res@info[[k]][[3]][2]
	 }
	 frob <- rbind(frobBC,frobSVD)
	 rownames(frob) <- c("Frobenius distance S4VD:","Frobenius distance SVD:")
	 colnames(frob) <- paste("BC", 1:res@Number)
	 print(frob)
	 invisible(frob)
}