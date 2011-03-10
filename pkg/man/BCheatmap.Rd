\name{BCheatmap}
\alias{OverlapHeat}
%\alias{s4vd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Overlap Heatmap for the visualization of overlapping biclusters
}
%\description{
%
%}
\usage{
BCheatmap(X,res,cexR=.75,cexC=.75,axisR=FALSE,axisC=TRUE,heatcols = diverge_hcl(25, h = c(260, 0), c = 80, l = c(30, 100), power = 1.5,gamma = 2.4, fixup = TRUE),
		clustercols= rainbow_hcl(res@Number, c = 100, l = 50),
		allrows=F,
		allcolumns=T
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
the data matrix  
}
  \item{res}{
the biclustering result
}
  \item{cexR}{
relativ font size of the row labels 
}
  \item{cexR}{
relativ font size of the column labels
}
  \item{axisR}{
if TRUE the row labels will be plotted
}
  \item{axisC}{
if TRUE the column labels will be plotted
}
  \item{heatcols}{
a character vector specifing the heatcolors 
}
  \item{axisC}{
a character vector specifing the colors of the rectangles that indicate the rows and columns that belong to a bicluster
}
  \item{allrows}{
if FALSE only the rows assigned to any bicluster will be plotted
}
  \item{allcolumns}{
if FALSE only the columns assigned to any bicluster will be plotted
}
}
%\details{
%plots the stability path the if savepath=TRUE. If the pointwise error control was used or if savepath=FALSE the final selection probabilities for the rows and the columns will be plotted.  
%}
%\value{
%}
%\references{
%
%kommt noch
%}
\author{
Martin Sill
}
%\note{
%kommt noch
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%kommt noch
%}
%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
