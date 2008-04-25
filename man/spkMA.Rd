\name{spkMA}
\alias{spkMA}
\title{MA Plots}
\description{
  Plots log-ratios (M) vs. average log expression (A) for a
  SpikeInExpressionSet object. 
}
\usage{
spkMA(object, spkSlopeOut, fc=2, tol=3, compare=NULL, label=NULL,
xaxis=NULL, ylim=NULL, outlier=1, reduce=TRUE)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{spkSlopeOut}{the output from the spkSlope function}
  \item{fc}{the fold change of interest}
  \item{tol}{the precision (number of digits after decimal point) in fc}
  \item{compare}{a matrix with two columns containing the arrays that
    should be compared when making fold change plots. This is used in
    particular for two-color arrays.}
  \item{label}{a character string to insert into the plot title}
  \item{xaxis}{comparisons to be made. Options are: Bg-Null LL,
    Bg-Null MM, Bg-Null HH, S-Null LL, S-Null MM, S-Null HH,
    LL, ML, MM, HL, HM, HH. If NULL, all comparisons will be made.}
  \item{ylim}{limits of y-axis}
  \item{outlier}{log fold change cut-off for outliers}
  \item{reduce}{if TRUE some points are removed from the background to
    speed plotting}
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
