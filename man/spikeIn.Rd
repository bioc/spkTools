\name{spikeIn}
\alias{spikeIn}
\title{Retrieve Nominal Concentrations from SpikeInExpressionSet objects}
\description{
  This function accesses the nominal concentrations of assay data stored
  in a SpikeInExpressionSet object.
}
\usage{
spikeIn(object)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
}
\value{
  returns a matrix of nominal concentration values.
}
\author{Matthew N. McCall}
\seealso{exprs}
\examples{
data(affy)
affySpikes <- spikeIn(affy)
}
\keyword{manip}

