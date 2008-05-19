\name{plotSpkBox}
\alias{plotSpkBox}
\title{Boxplots of Fold Changes Calculated by spkBox}
\description{
  Plots boxplots of the data resulting from a call to spkBox.
}
\usage{
plotSpkBox(boxs, fc=2, box.names=NULL, ...)
}
\arguments{
  \item{boxs}{the output of a call to spkBox}
  \item{fc}{expected fold change}
  \item{box.names}{names to be printed below each boxplot}
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
plotSpkBox(affyBox) 
}
\keyword{hplot}
