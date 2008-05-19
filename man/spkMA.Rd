\name{spkMA}
\alias{spkMA}
\title{MA Plots}
\description{
  Plots log-ratios (M) vs. average log expression (A) for a
  SpikeInExpressionSet object. 
}
\usage{
spkMA(object, spkSlopeOut, fc=2, tol=3, label=NULL, ylim=NULL,
outlier=1, reduce=TRUE, plot.legend=TRUE)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{spkSlopeOut}{the output from the spkSlope function}
  \item{fc}{the fold change of interest}
  \item{tol}{the precision (number of digits after decimal point) in fc}
  \item{label}{a character string to insert into the plot title}
  \item{ylim}{limits of y-axis}
  \item{outlier}{log fold change cut-off for outliers}
  \item{reduce}{if TRUE some points are removed from the background to
    speed plotting}
  \item{plot.legend}{if TRUE a legend is plotted}
}
\value{
  The MA plot is produced.
}
\author{Matthew N. McCall}
\examples{
data(affy)
affySlope <- spkSlope(affy)
spkMA(affy, affySlope)
}
\keyword{hplot}
