\name{BCs4vd}
\alias{BCs4vd}
\alias{s4vd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Biclustering via sparse singular value decomposition incorporating stability selection
}
\description{
The function performs the a biclustering of the data matrix by sparse singular value decomposition with a nested stability selection.
}
\usage{
\S4method{biclust}{matrix,BCs4vd}(x, method=BCs4vd(),
		steps = 100,
		pcerv = 0.05,
		pceru = 0.05,
		ss.thr = c(0.6,0.65),
		size = 0.632,
		gamm = 0,
		iter = 100,
		nbiclust = 10,
		merr = 0.0001,
		cols.nc=FALSE,
		rows.nc=TRUE,
		row.overlap=TRUE,
		col.overlap=TRUE,
		row.min=4,
		col.min=4,
		pointwise=TRUE,
		start.iter=0,
		savepath=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The matrix to be clustered.
}
  \item{steps}{
Number of subsamples used to perform the stability selection.
}
  \item{pcerv}{
Per comparsion wise error rate to control the number of falsely selected right singular vector coefficients (columns/samples).
}
  \item{pceru}{
Per comparsion wise error rate to control the number of falsely selected left singular vector coefficients (rows/genes).
}
  \item{ss.thr}{
Range of the cutoff threshold (relative selection frequency) for the stability selection.
}
  \item{size}{
Size of the subsamples used to perform the stability selection.  
}
  \item{gamm}{
Weight parameter for the adaptive LASSO, nonnegative constant (default = 0, LASSO).
}
  \item{iter}{
Maximal number of iterations to fit a single bicluster.
}
  \item{nbiclust}{
Maximal number of biclusters. 
}
  \item{merr}{
Threshold to decide convergence. 
}
  \item{cols.nc}{
Allow for negative correlation of columns (samples) over rows (genes).
}
  \item{rows.nc}{
Allow for negative correlation of rows (genes) over columns (samples).
}
  \item{row.overlap}{
Allow rows to overlap between biclusters. 
}
  \item{col.overlap}{
Allow columns to overlap between biclusters. 
}
  \item{row.min}{
Minimal number of rows.
}
  \item{col.min}{
Minimal number of columns.
}
  \item{pointwise}{
If TRUE performs a fast pointwise stability selection instead of calculating the complete stability path.  
}
  \item{start.iter}{
Number of starting iterations in which the algorithm is not allowed to converge. 
}
  \item{savepath}{
Saves the stability path in order plot the path with the stabpathplot function.
Note that pointwise needs to be TRUE to save the path. For extreme high dimensional data sets (e.g. the lung cancer example) the resulting
biclust object may exceed the available memory.
}
}
\value{
Returns an object of class \code{Biclust}.
}
\author{
Martin Sill
\email{m.sill@dkfz.de}
}

\seealso{ \code{\link{biclust}}, \code{\link{Biclust}}}

\examples{

}

\keyword{cluster}
\keyword{classif}
\keyword{svd}
\keyword{stability selection}
