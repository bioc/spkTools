\name{spkAll}
\alias{spkAll}
\title{Spike-in Functions Wrapper}
\description{
  A wrapper for the functions contained in the spkTools package, which
  calls each function.
}
\usage{
spkAll(object, label, model=expr~spike+probe+array, fc=NULL, tol=3,
xrngs=NULL, yrngs=NULL, cuts=c(.6,.99), potQuantile=.995,
gnn=c(25,100,10000), pch=".", output="eps")
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{label}{a character string to insert into the graphs and tables
    produced}
  \item{model}{model to be passed to spkAnova}
  \item{fc}{the fold change for which fold change plots will be produced}
  \item{tol}{the number of digits after the decimal point in fc}
  \item{xrngs}{ranges for the x-axis of each plot. d=density, s=slope,
    v=box, m=M vs A}
  \item{yrngs}{ranges for the y-axis of each plot. d=density, s=slope,
    v=box, m=M vs A}
  \item{cuts}{quantiles used to make the low, medium, and high bins}
  \item{potQuantile}{the desired quantile to compute the probability of
    being above}
  \item{gnn}{a vector of 3 numbers passed to spkGNN: the desired number
    of true positives, the number of truly expressed genes, and the
    number of truly unexpressed genes}
  \item{pch}{plotting point to be used in spkSlope}
  \item{output}{the format in which to save the plots produced. Options
    are "pdf" and "eps"}
}
\value{
  The full complement of plots and tables described in the vignette are
  created and saved in the current working directory.
}
\author{Matthew N. McCall}
\examples{
data(affy)
spkAll(affy, label="affy", fc=2)
}
\keyword{hplot}
\keyword{manip}


