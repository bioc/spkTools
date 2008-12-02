\name{SpikeInExpressionSet-class}
\docType{class}
\alias{class:SpikeInExpressionSet}
\alias{SpikeInExpressionSet-class}
\alias{SpikeInExpressionSet}
\alias{SpikeInExpressionSet-methods}
\alias{initialize,SpikeInExpressionSet-method}
\alias{spikeIn,SpikeInExpressionSet-method}
\alias{spikeIn<-,SpikeInExpressionSet,matrix-method}
\alias{spkSplit,SpikeInExpressionSet-method}
\title{
  Class to Contain and Describe High-Throughput Expression Level Assays
  with Spike-in Data}
\description{This is a class representation for spike-in expression
  data. \code{SpikeInExpressionSet} class is dervied from
  \code{\link{ExpressionSet}}, and requires a matrix names \code{exprs}
  and a matrix named \code{spikeIn}.}
\section{Extends}{
  Extends class \code{\link{ExpressionSet}}.
  }
\section{Creating Objects}{

  \code{new("SpikeInExpressionSet")}

  \code{new("SpikeInExpressionSet",
    phenoData = new("AnnotatedDataFrame"),
    featureData = new("AnnotatedDataFrame"),
    experimentData = new("MIAME"),
    annotation = character(0),
    exprs = new("matrix"),
    spikeIn = new("matrix"))
  }  

  This creates a \code{SpikeInExpressionSet} with \code{assayData}
  implicitly created to contain \code{exprs} and
  \code{spikeIn}. Additional named matrix arguements with the same
  dimensions as \code{exprs} are added to \code{assayData}; the row
  and column names of these additional matricies should match those of
  \code{exprs} and \code{spikeIn}.

  \code{new("SpikeInExpressionSet",
    assayData = assayDataNew(exprs=new("matrix"),spikeIn=new("matrix")),  
    phenoData = new("AnnotatedDataFrame"),
    featureData = new("AnnotatedDataFrame"),
    experimentData = new("MIAME"),
    annotation = character(0),
  }

  This creates a \code{SpikeInExpressionSet} with \code{assayData}
  provided explicitly. In this form, the only required named argument is
  \code{assayData}.
}
    
\section{Slots}{
  Inherited from \code{ExpressionSet}:
  \describe{
    \item{\code{assayData}:}{Contains matrices with equal dimensions,
      and with column number equal to
      \code{nrow(phenoData)}. \code{assayData} must contain a matrix
      \code{exprs} and a matrix \code{spikeIn} with rows representing
      features and columns representing samples.}  
    \item{\code{phenoData}:}{See \code{\link{eSet}}}
    \item{\code{annotation}}{See \code{\link{eSet}}}
    \item{\code{featureData}}{See \code{\link{eSet}}}
    \item{\code{experimentData}:}{See \code{\link{eSet}}} 
  }
}

\section{Methods}{
  Class-specific methods:
  \describe{
    \item{\code{spikeIn(SpikeInExpressionSet)},
      \code{spikeIn(SpikeInExpressionSet)<-}}{Access and set elements
      named \code{spikeIn} in the \code{AssayData-class} slot.}
    \item{\code{spkSplit(SpikeInExpressionSet)}}{creates two SpikeInExpressionSet objects
      -- one with the spike-in probes and one with the non-spike-in
      probes.}
  }

  For derived methods (see \code{\link{ExpressionSet}}).
}

\seealso{
  \code{\link{eSet-class}}, \code{\link{ExpressionSet-class}}. 
} 

\examples{
# create an instance of SpikeInExpressionSet
new("SpikeInExpressionSet")

new("SpikeInExpressionSet", exprs=matrix(runif(1000), nrow=100), spikeIn=matrix(rep(1:10,100), nrow=100))

# class specific methods
data(affy)
affySpikes <- spikeIn(affy)
affySplit <- spkSplit(affy)
}

\keyword{classes}
