\name{spkSlope}
\alias{spkSlope}
\title{Signal Detect Slope Plot}
\description{
  Plots observed expression vs. nominal concentration. The overall
  regression slope, as well as, regression slopes for low, medium, and
  high bins are computed and the regression lines plotted. 
}
\usage{
spkSlope(object, label = NULL, cuts=c(.6,.99), ...)
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{label}{a character string to insert into the plot title}
  \item{cuts}{quantiles used to make the low, medium, and high bins}
  \item{...}{arguments passed to the plot function}
}
\details{
  The bins are created by computing the proportion of non-spike-in genes
  with expression values less than or equal to the average expression
  value at each nominal concentration. Using the default value of cuts,
  the high bin contains nominal concentrations with 99 precent or more
  of the non-spike-in expression 
  values lower than it. The medium bin contains nominal concentrations
  with between 60 and 99 precent of the non-spike-in expression values
  lower than it. The low bin contains nominal concentrations with less
  than 60 percent of the non-spike-in expression values lower than it.
}
\value{
  \item{avgExp}{average expression at each nominal concentration}
  \item{slopes}{the regression slopes - overall and for each bin}
  \item{breaks}{which spike-in levels fall in each bin}
  \item{brkpts}{the expression value of the cut points between bins}
  \item{prop}{the proportion of non-spike-in probes with expression less
  than the average expression at each nominal concentration}
}
\author{Matthew N. McCall}
\examples{
data(affy)
spkSlope(affy)
}
\keyword{hplot}
\keyword{manip}

