\name{spkVar}
\alias{spkVar}
\title{Spike-in Variance}
\description{
  Compute an estimate of the standard deviation in expression at each nominal concentration. 
}
\usage{
spkVar(object)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
}
\value{
  a matrix containing spike-in levels and corresponding MADs.
}
\author{Matthew N. McCall}
\examples{
data(affy)
spkVar(affy)
}
\keyword{manip}

