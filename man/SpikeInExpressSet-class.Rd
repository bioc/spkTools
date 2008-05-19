\name{SpikeInExpressionSet-class}
\docType{class}
\alias{SpikeInExpressionSet-class}
\alias{SpikeInExpressionSet,ANY}
\alias{SpikeInExpressionSet}
\alias{initialize,SpikeInExpressionSet-class}
\alias{spikeIn,SpikeInExpressionSet-class}
\alias{spkSplit,SpikeInExpressionSet-class}
\title{SpikeInExpressionSet-class}
\description{This is a class representation for spike-in expression
  data. The main components are the intensities and nominal
  concentrations from multiple arrays. It extends
  \code{\link[Biobase:class.eSet]{eSet}}.} 
\section{Objects from the Class}{
  Objects can be created by creating a new object of type SpikeInExpressionSet.
}
\section{Slots}{
  \describe{
    \item{\code{assayData}:}{Object of class \code{AssayData}
      containing the raw data, which will be at minimum a matrix of
      intensity values (exprs) and a matrix of nominal concentrations
      (spikeIn).}  
    \item{\code{phenoData}:}{Object of class \code{AnnotatedDataFrame}
      containing phenotypic data for the samples.}
    \item{\code{annotation}}{A character string identifying the
      annotation used for the \code{ExpressionSet} instance.}
    \item{\code{featureData}}{Object of class \code{AnnotatedDataFrame}
      containing feature-level (e.g., probeset-level) information.}
    \item{\code{experimentData}:}{Object of class "MIAME" containing
      experiment-level information.} 
  }
}
\section{Extends}{
Class \code{"ExpressionSet"}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(object = "SpikeInExpressionSet",
	method="character")}} used to create a new SpikeInExpressionSet
    object.
    \item{spikeIn}{\code{signature(object = "SpikeInExpressionSet",
	method="character")}} extract the nominal concentration matrix
    from a SpikeInExpressionSet object.
    \item{spkSplit}{\code{signature(object = "SpikeInExpressionSet",
	method="character")}} creates two SpikeInExpressionSet objects
    -- one with the spike-in probes and one with the non-spike-in
    probes. 
  }
}
\note{This class is better described in the vignette.}
\seealso{\code{\link[Biobase:class.eSet]{eSet}}} 

\examples{
data(affy)
affySpikes <- spikeIn(affy)
affySplit <- spkSplit(affy)
}

\keyword{classes}
