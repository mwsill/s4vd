setClass('BCs4vd',
		contains = 'BiclustMethod',
		prototype = prototype(
				biclustFunction = function(
						X,
						steps = 100,
						pcerv = 0.1,
						pceru = 0.1,
						ss.thr = c(0.6,0.65),
						size = 0.5,
						gamm = 0,
						iter = 100,
						nbiclust = 10,
						merr = 10^(-3),
						cols.nc=TRUE,
						rows.nc=TRUE,
						row.overlap=TRUE,
						col.overlap=TRUE,
						row.min=1,
						col.min=1,
						pointwise=TRUE,
						start.iter=3,
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

