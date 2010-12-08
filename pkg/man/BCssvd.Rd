\name{BCssvd}
\alias{BCssvd}
\alias{ssvd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Biclustering via sparse singular value decomposition 
}
\description{
The function performs a biclustering of the data matrix by sparse singular value decomposition.
\usage{
ssvd(X,threu = 1, threv = 1, gamu = 0, gamv =0,  u0 = svd(X)$u[,1], v0 = svd(X)$v[,1], merr = 10^(-4), niter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
the matrix to be clustered
}
  \item{threu}{
type of penalty (thresholding rule) for the left singular vector,
				 1 = (Adaptive) LASSO (default)
                 2 = hard thresholding
}
  \item{threv}{
type of penalty (thresholding rule) for the right singular vector,
                 1 = (Adaptive) LASSO (default)
                 2 = hard thresholding	
}
  \item{gamu}{
weight parameter in Adaptive LASSO for the left singular vector, nonnegative constant (default = 0, LASSO)
}
  \item{gamv}{
weight parameter in Adaptive LASSO for the right singular vector, nonnegative constant (default = 0, LASSO)
}
  \item{u0}{
initial left singular vector                
}
  \item{v0}{
initial right singular vector   
}
  \item{merr}{
threshold to decide convergence
}
  \item{niter}{
maximum number of iterations
}
}

\value{
Returns an Biclust object.
}
\references{
Mihee Lee, Haipeng Shen, Jianhua Z. Huang and J. S. Marron1
"Biclustering via Sparse Singular Value Decomposition", Biometrics, 2010 
}
\author{
Adaptation of original code from Mihee Lee by Martin Sill
\email{m.sill@dkfz.de}
}

\seealso{ \code{\link{biclust}}, \code{\link{Biclust}}}

\examples{

}

\keyword{cluster}
\keyword{classif}
\keyword{svd}
