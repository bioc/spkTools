\name{spkSplit}
\alias{spkSplit}
\title{Split a SpikeInExpressionSet Object}
\description{
  Takes a SpikeInExpressionSet object and creates two new
  SpikeInExpressionSet objects -- one with only the spike-in probes and
  one with the non-spike-in probes.
}
\usage{
spkSplit(object)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
}
\value{
  A list with the following two elements:
  \item{s}{SpikeInExpressionSet object of spike-in probes}
  \item{ns}{SpikeInExpressionSet object of non-spike-in probes}
}
\author{Matthew N. McCall}
\examples{
data(affy)
affySplit <- spkSplit(affy)
}
\keyword{manip}

