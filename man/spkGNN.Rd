\name{spkGNN}
\alias{spkGNN}
\title{Genes Needed to Detect N True Positives}
\description{
  Computes the number of genes one would need to consider to obtain a
  given number of truly positive genes if one considered genes in order
  of decreasing observed fold change.
}
\usage{
spkGNN(n, n.expr, n.unexpr, AccuracySlope, AccuracySD, nullfc)
}
\arguments{
  \item{n}{the desired number of true positives}
  \item{n.expr}{the actual number of truly expressed genes}
  \item{n.unexpr}{the actual number of truly unexpressed genes}
  \item{AccuracySlope}{the signal detect slope from the spkSlope function}
  \item{AccuracySD}{the standard deviation of the signal detect slope
    from the spkAccSD function}
  \item{nullfc}{a vector of null fold changes from the spkBox function}  
}
\value{
  This function returns the expected number of genes one would have to
  consider to obtain N true positives under the given conditions.
}
\author{Matthew N. McCall}
\examples{
data(affy)
spkSlopeOut <- spkSlope(affy)
spkBoxOut <- spkBox(affy, spkSlopeOut, fc=2)
AccuracySlope <- round(spkSlopeOut$slope[-1], digits=2)
AccuracySD <- round(spkAccSD(affy, spkSlopeOut), digits=2)
spkGNN(n=25, n.expr=100, n.unexpr=10000, AccuracySlope[2],
AccuracySD[2], spkBoxOut[[2]])
}
\keyword{manip}


