\name{spkAccSD}
\alias{spkAccSD}
\title{Accuracy Standard Deviation}
\description{
  Estimates the standard deviation for spike-ins at the lowest possible
  fold change in each bin.
}
\usage{
spkAccSD(object, spkSlopeOut, tol=3)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{spkSlopeOut}{the output from the spkSlope function}
  \item{tol}{number of digits after decimal point}
}
\value{
  returns the median absolute deviation (MAD) for each bin. 
}
\author{Matthew N. McCall}
\examples{
data(affy)
affySlope <- spkSlope(affy)
spkAccSD <- spkAccSD(affy, affySlope)
}
\keyword{manip}
