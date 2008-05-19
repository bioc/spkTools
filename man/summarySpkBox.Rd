\name{summarySpkBox}
\alias{summarySpkBox}
\title{Summary of Fold Changes Calculated by spkBox}
\description{
  Prints a summary table of the data resulting from a call to spkBox.
}
\usage{
summarySpkBox(boxs)
}
\arguments{
  \item{boxs}{the output of a call to spkBox}
}
\value{
  A dataframe with 2 columns: the mean fold change and the median
  average distance of the fold changes.
}
\author{Matthew N. McCall}
\examples{
data(affy)
affySlope <- spkSlope(affy)
affyBox <- spkBox(affy, affySlope)
plotSpkBox(affyBox) 
}
\keyword{manip}
