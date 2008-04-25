\name{spkAll}
\alias{spkAll}
\alias{spkFuncs}
\title{Spike-in Functions Wrapper}
\description{
  A wrapper for the functions contained in the spkTools package, which
  calls each function for each sample.
}
\usage{
spkAll(object, label, model=expr~spike+probe+array, fc=NULL, tol=3,
compare=NULL, xrngs=NULL, yrngs=NULL, grps=NULL, cuts=c(.6,.99),
potQuantile=.995, xaxis=NULL, pch=".", output="eps")
spkFuncs(object, label, model=expr~spike+probe+array, fc=NULL, tol=3,
compare=NULL, xrngs=NULL, yrngs=NULL, grps=NULL, cuts=c(.6,.99),
potQuantiles=.995, xaxis=NULL, pch=".", output="eps")
}
\arguments{
  \item{object}{a SpikeInExpressionSet object}
  \item{label}{a character string to insert into the graphs and tables
    produced}
  \item{model}{model to be passed to spkAnova}
  \item{fc}{the fold change for which fold change plots will be produced}
  \item{tol}{the number of digits after the decimal point in fc}
  \item{compare}{a matrix with two columns containing the arrays that
    should be compared when making fold change plots. This is used in
    particular for two-color arrays.}
  \item{xrngs}{ranges for the x-axis of each plot. d=density, s=slope,
    v=box, m=M vs A}
  \item{yrngs}{ranges for the y-axis of each plot. d=density, s=slope,
    v=box, m=M vs A}
  \item{grps}{comparisons to be made in spkBox. Options are: Bg-Null LL,
    Bg-Null MM, Bg-Null HH, S-Null LL, S-Null MM, S-Null HH,
    LL, ML, MM, HL, HM, HH.}
  \item{cuts}{quantiles used to make the low, medium, and high bins}
  \item{potQuantile}{the desired quantile to compute the probability of
    being above}
  \item{xaxis}{comparisons to be made in spkMA. Options are: Bg-Null LL,
    Bg-Null MM, Bg-Null HH, S-Null LL, S-Null MM, S-Null HH,
    LL, ML, MM, HL, HM, HH. If NULL, all comparisons will be made.}
  \item{pch}{plotting point to be used in spkSlope}
  \item{output}{the format in which to save the plots produced. Options
    are "pdf" and "eps".}
}
\details{
  If there are multiple samples in the SpikeInExpressionSet object that
  can be accessed via pData(object)\$sample, the plots and tables will be
  computed for each sample separately. 
}
\value{
  A number of plots and tables are created and saved in the current
  working directory.
}
\author{Matthew N. McCall}
\examples{
data(affy)
spkAll(affy, label="affy", fc=2)
}
\keyword{hplot}


