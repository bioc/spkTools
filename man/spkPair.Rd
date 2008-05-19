\name{spkPair}
\alias{spkPair}
\title{Pairwise Comparisons for Spike-in Genes}
\description{
  Compute log-ratios among spike-in genes. 
}
\usage{
spkPair(object)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
}
\value{
  An array containing either log-ratios (M), average log expression
  (A), and nominal concentrations (N1 & N2). Dimension one is genes,
  dimension two is array pairings, dimension three is M, A, N1, and N2. 
}
\author{Matthew N. McCall}
\examples{
data(affy)
affyPair <- spkPair(affy)
}
\keyword{manip}
