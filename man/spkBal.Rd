\name{spkBal}
\alias{spkBal}
\title{Quantify Microarry Spike-in Design Imbalance}
\description{
  Computes the imbalance of a microarray spike-in design due to probes
  and arrays.
}
\usage{
spkBal(object)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
}
\value{
  The probe and array imbalances.
}
\references{Wu, Chien-Fu, Iterative Construction of Nearly Balanced
  Assignments I: Categorical Covariates. Technometrics, Vol. 23,
  No. 1. (Feb, 1981), pp. 37-44.}
\author{Matthew N. McCall}
\examples{
data(affy)
spkBal(affy)
}
\keyword{manip}

