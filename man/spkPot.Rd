\name{spkPot}
\alias{spkPot}
\title{Probability of being in the Top}
\description{
  Compute the probability that a spike-in with a nominal fold change of
  2 appears in the the top 0.5\% (default) of log-ratios. 
}
\usage{
spkPot(object, spkSlopeOut, sig, SD, precisionQuantile)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{spkSlopeOut}{the output from the spkSlope function}
  \item{sig}{the signal detect slopes from a call to spkSlope}
  \item{SD}{the standard deviation from spkAccSD}
  \item{precisionQuantile}{the desired quantile to compute the
    probability of being above}
}
\value{
  A vector of probabilities for each ALE strata. 
}
\author{Matthew N. McCall}
\examples{
data(affy)
affySlope <- spkSlope(affy)
affyAccSD <- spkAccSD(affy, affySlope)
spkPot(affy, affySlope, affySlope$slopes, affyAccSD, .995)
}
\keyword{manip}
