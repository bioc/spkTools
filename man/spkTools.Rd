\name{spkTools}
\alias{spkTools-package}
\alias{spkTools}
\docType{package}
\title{
Tools for Spike-in Data Analysis and Visualization
}
\description{
  A collection of functions to examine microarray datasets that
  include spike-ins. In particular, it allows one to explore the
  distribution of spike-ins within the range of possible expression
  values, the relationship between nominal concentration and expression,
  and the relationship between expected and observed fold change for
  different levels of comparison.
}
\details{
\tabular{ll}{
Package: \tab spkTools\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2007-10-9\cr
License: \tab GPL version 2 or newer\cr
}

}
\author{
Matthew N. McCall

Maintainer: Matthew N. McCall <mmccall@jhsph.edu>
}
\keyword{ package }
\examples{
## The Three Plots
data(affy)
par(mfrow=c(2,2))
affySlope <- spkSlope(affy)
spkDensity(affy, affySlope)
spkBox(affy, affySlope)

## The Full Wrapper
data(affy)
spkAll(affy, label="Affymetrix", fc=2)
}
