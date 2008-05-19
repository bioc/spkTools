\name{spkPairNS}
\alias{spkPairNS}
\title{Pairwise Comparisons for Non-Spike-in Genes}
\description{
  Compute log-ratios among non-spike-in genes. 
}
\usage{
spkPairNS(object, output="M")
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{output}{what to return; either "M" for log-ratios or "A" for
    average log expression.}
}
\value{
  A matrix containing either log-ratios (M) or average log expression
  (A). Rows are genes and columns are array pairings.
}
\author{Matthew N. McCall}
\examples{
data(affy)
affyPairNS <- spkPairNS(affy)
}
\keyword{manip}
