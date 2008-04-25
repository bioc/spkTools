\name{spkPairNS}
\alias{spkPairNS}
\title{Pairwise Comparisons for Non-Spike-in Genes}
\description{
  Compute log-ratios among non-spike-in genes. 
}
\usage{
spkPairNS(object, compare=NULL, output="M")
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{compare}{a matrix with two columns of which arrays to
    compare. If NULL, all possible combinations are made.}
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
