\name{spkBox}
\alias{spkBox}
\title{Fold Change Calculations}
\description{
  A function to calculate the log-ratios stratified by which ALE groups
  yield the comparison. They are stratified by which bins are being
  compared to produce the given fold change.
}
\usage{
spkBox(object, spkSlopeOut, fc = 2, tol = 3, compare = NULL,
grps=NULL)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{spkSlopeOut}{the output of the spkSlope function}
  \item{fc}{the fold change of interest}
  \item{tol}{the precision (number of digits after decimal point) in fc}
  \item{compare}{a matrix with two columns containing the arrays that
    should be compared when making fold change plots. This is used in
    particular for two-color arrays.}
  \item{grps}{comparisons to be made. Options are: Bg-Null LL,
    Bg-Null MM, Bg-Null HH, S-Null LL, S-Null MM, S-Null HH,
    LL, ML, MM, HL, HM, HH. If NULL, all comparisons will be made.}
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

