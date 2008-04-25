\name{spkPair}
\alias{spkPair}
\title{Pairwise Comparisons for Spike-in Genes}
\description{
  Compute log-ratios among spike-in genes. 
}
\usage{
spkPairNS(object, compare=NULL)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{compare}{a matrix with two columns of which arrays to
    compare. If NULL, all possible combinations are made.}
}
\value{
  An array containing either log-ratios (M), average log expression
  (A), and nominal concentrations (N1 & N2). Dimension one is genes,
  dimension two is array pairings, dimension three is M, A, N1, and N2. 
}
\author{Matthew N. McCall}
\examples{
data(affy)
affyPair <- spkPairNS(affy)
}
\keyword{manip}
