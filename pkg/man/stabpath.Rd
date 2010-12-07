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
s4vd(X, r.steps = 500, c.steps = 500, r.err = 20, c.err = 10, iter = 50, ss.thr = 0.8, size = 0.5, weak = 0.2, r.overlap = TRUE, c.overlap = FALSE, r.negcorr = TRUE, c.negcorr = FALSE, nbiclust = 50, merr = 0.01, mc.cores = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
the matrix to be clustered
}
  \item{r.steps}{
number of subsamples to calculate the relative selection frequencies of the left singular vector coefficients that correspond to the rows of X
}
  \item{c.steps}{
number of subsamples to calculate the relative selection frequencies of the right singular vector coefficients that correspond to the columns of X	
}
  \item{r.err}{
expected number of falsely selected left singular vector coefficients 	
}
  \item{c.err}{
expected number of falsely selected right singular vector coefficients
}
  \item{iter}{
maximal number of iterations to fit a single bicluster
}
  \item{ss.thr}{
cutoff threshold (relative selection frequency) for the stability selection
}
  \item{size}{
size of the subsamples for the stability selection
}
  \item{weak}{
weakness parameter for the randomised lasso
}
  \item{r.overlap}{
allow bicluster to be row overlapping
}
  \item{c.overlap}{
allow bicluster to be column overlapping
}
  \item{r.negcorr}{
allow for negative correlation of rows (genes) over columns (samples)
}
  \item{c.negcorr}{
allow for negative correlation of columns (samples) over rows (genes)
}
  \item{nbiclust}{
maximal number of biclusters
}
  \item{merr}{
convergence parameter
}
  \item{mc.cores}{
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
