\name{spkBox}
\alias{spkBox}
\title{Fold Change Calculations}
\description{
  A function to calculate the log-ratios stratified by which ALE groups
  yield the comparison. They are stratified by which bins are being
  compared to produce the given fold change.
}
\usage{
spkBox(object, spkSlopeOut, fc = 2, tol = 3, reduce=TRUE)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{spkSlopeOut}{the output of the spkSlope function}
  \item{fc}{the fold change of interest}
  \item{tol}{the precision (number of digits after decimal point) in fc}
  \item{reduce}{if TRUE the number of points plotted in the null bins is reduced}
}
\details{
  This function requires the output of spkSlope.
}
\value{
  A list with the log-ratios separated by ALE strata comparison.
}
\author{Matthew N. McCall}
\examples{
data(affy)
affySlope <- spkSlope(affy)
spkBox(affy,affySlope)
}
\keyword{manip}

