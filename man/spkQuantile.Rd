\name{spkQuantile}
\alias{spkQuantile}
\title{Empirical Quantiles}
\description{
  An internal function called by spkSlope.
}
\usage{
spkQuantile(amt, avgE, ens, p)
}
\arguments{
  \item{amt}{a vector of nominal concentrations}
  \item{avgE}{the observed average expression corresponding to each
    nominal concentration}
  \item{ens}{the average expression across arrays of unexpressed genes}
  \item{p}{the quantiles to make the bins}
}
\author{Matthew N. McCall}
\examples{
data(affy)
affySlope <- spkSlope(affy)
}
\keyword{manip}
