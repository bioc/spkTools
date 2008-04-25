\name{plot.spkBox}
\alias{plot.spkBox}
\title{Boxplots of Fold Changes Calculated by spkBox}
\description{
  Plots boxplots of the data resulting from a call to spkBox.
}
\usage{
plot.spkBox(boxs, fc=2, ...)
}
\arguments{
  \item{boxs}{the output of a call to spkBox}
  \item{fc}{expected fold change}
  \item{...}{parameters passed to boxplot}
}
\value{
  Boxplots are produced.
}
\author{Matthew N. McCall}
\examples{
data(affy)
affySlope <- spkSlope(affy)
affyBox <- spkBox(affy, affySlope)
plot.spkBox(affyBox) 
}
\keyword{hplot}
