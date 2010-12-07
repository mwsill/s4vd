\name{BCs4vd}
\alias{BCs4vd}
\alias{s4vd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Biclustering via sparse singular value decomposition incorporating stability selection
}
\description{
The function performs the a biclustering of the data matrix by sparse singular value decomposition with a nested stability selection
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
the matrix to be clustered
}
  \item{steps}{
number of subsamples used to perform the stability selection
}
  \item{pcerv}{
per comparsion wise error rate to control the number of falsely selected right singular vector coefficients (columns/samples)
}
  \item{pceru}{
per comparsion wise error rate to control the number of falsely selected left singular vector coefficients (rows/genes)
}
  \item{ss.thr}{
range of the cutoff threshold (relative selection frequency) for the stability selection
}
  \item{size}{
size of the subsamples used to perform the stability selection  
}
  \item{gamm}{
weight parameter for the adaptive LASSO, nonnegative constant (default = 0, LASSO)
}
  \item{iter}{
maximal number of iterations to fit a single bicluster
}
  \item{nbiclust}{
maximal number of biclusters 
}
  \item{merr}{
threshold to decide convergence 
}
  \item{cols.nc}{
allow for negative correlation of columns (samples) over rows (genes)
}
  \item{rows.nc}{
allow for negative correlation of rows (genes) over columns (samples)
}
  \item{row.overlap}{
allow rows to overlap between biclusters 
}
  \item{col.overlap}{
allow columns to overlap between biclusters 
}
  \item{row.min}{
minimal number of rows
}
  \item{col.min}{
minimal number of columns
}
  \item{pointwise}{

}
  \item{start.iter}{
number 
}
  \item{savepath}{
number of cores for parallelization
}
}
\details{
kommt noch
}
\value{
kommt noch
}
\references{
kommt noch
}
\author{
kommt noch
}
\note{
kommt noch
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
kommt noch
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
data(BicatYeast)
X <- scale(BicatYeast)
set.seed(12345)
res <- biclust(X,method=BCs4vd,r.err=50,c.err=5,r.steps=100,c.steps=100,nbiclust=5)
## The function is currently defined as

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
