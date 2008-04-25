\name{spkDensity}
\alias{spkDensity}
\title{Spike-in Density Plot}
\description{
  A density plot of the non-spike-in expression with a rug of the average expression at each spike-in level.
}
\usage{
spkDensity(object, spkSlopeOut, cuts=TRUE, label = NULL, ...)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{spkSlopeOut}{the output from the spkSlope function}
  \item{cuts}{if TRUE vertical lines are drawn at the expression values
    separating low vs medium and medium vs high ALE strata}
  \item{label}{a character string to insert into the plot title}
  \item{...}{arguments passed to the plot function}
}
\details{
  This function requires the output of spkSlope.
}
\value{
  Density plot is produced.
}
\author{Matthew N. McCall}
\examples{
data(affy)
affySlope <- spkSlope(affy)
spkDensity(affy,affySlope)
}
\keyword{hplot}

