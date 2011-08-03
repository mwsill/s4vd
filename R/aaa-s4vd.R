setClass('BCs4vd',
		contains = 'BiclustMethod',
		prototype = prototype(
				biclustFunction = function(
						X,
						steps = 100,
						pcerv = 0.05,
						pceru = 0.05,
						ss.thr = c(0.6,0.65),
						size = 0.632,
						gamm = 0,
						iter = 100,
						nbiclust = 10,
						merr = 10^(-4),
						cols.nc=FALSE,
						rows.nc=TRUE,
						row.overlap=TRUE,
						col.overlap=TRUE,
						row.min=4,
						col.min=4,
						pointwise=TRUE,
						start.iter=0,
						savepath=FALSE)
					{s4vd(X,steps,pcerv,pceru,ss.thr,size,gamm,
						  iter,nbiclust,merr,cols.nc,rows.nc,
								row.overlap,
								col.overlap,
								row.min,
								col.min,
								pointwise,
								start.iter,
								savepath)}))



BCs4vd <- function() {
	return(new('BCs4vd'))
}


setClass('BCssvd',
		contains = 'BiclustMethod',
		prototype = prototype(
				biclustFunction = function(X,...){ssvdBC(X,...)}))

BCssvd <- function() {
	return(new('BCssvd'))
}

