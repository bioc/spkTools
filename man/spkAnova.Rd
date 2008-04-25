\name{spkAnova}
\alias{spkAnova}
\title{Anova Model for Microarray Spike-in Data}
\description{
  Computes the mean squared errors of a microarray spike-in design due
  to concentration, probe, array, and error. 
}
\usage{
spkAnova(object, model=expr~spike+probe+array)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{model}{the anova model}
}
\value{
  A vector of the mean squared errors from the anova model.

}
\author{Matthew N. McCall}
\examples{
data(affy)
spkAnova(affy)
}
\keyword{manip}

