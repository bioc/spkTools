## generic functions
setGeneric("spikeIn", function(object) standardGeneric("spikeIn"))
setGeneric("spkSplit", function(object) standardGeneric("spkSplit"))
setGeneric("spkDensity", function(object, spkSlopeOut, cuts=TRUE, label=NULL, ...)
           standardGeneric("spkDensity"))
setGeneric("spkSlope", function(object, label=NULL,
                                cuts=c(.6,.99), ...)
           standardGeneric("spkSlope"))
setGeneric("spkPair",function(object, compare=NULL) standardGeneric("spkPair"))
setGeneric("spkPairNS",function(object, compare=NULL, output="M")
           standardGeneric("spkPairNS"))
setGeneric("spkVar", function(object) standardGeneric("spkVar"))
setGeneric("spkBox",function(object, spkSlopeOut, fc=2, tol=3, compare=NULL, grps=NULL)
           standardGeneric("spkBox"))
setGeneric("spkMA",function(object, spkSlopeOut, fc=2, tol=3, compare=NULL,
                            label=NULL, xaxis=NULL, ylim=NULL, outlier=1,
                            reduce=TRUE)
           standardGeneric("spkMA"))
setGeneric("spkFuncs", function(object, label, model=expr~spike+probe+array,
                                fc=NULL, tol=3, compare=NULL, xrngs=NULL,
                                yrngs=NULL, grps=NULL, cuts=c(.6,.99),
                                potQuantile=.995, xaxis=NULL,
                                pch=".", output="eps")
           standardGeneric("spkFuncs"))
setGeneric("spkAll", function(object, label, model=expr~spike+probe+array,
                              fc=NULL, tol=3, compare=NULL, xrngs=NULL,
                              yrngs=NULL, grps=NULL,
                              cuts=c(.6,.99), potQuantile=.995, xaxis=NULL,
                              pch=".", output="eps")
           standardGeneric("spkAll"))
setGeneric("spkBal",function(object) standardGeneric("spkBal"))
setGeneric("spkAnova", function(object, model=expr~spike+probe+array) standardGeneric("spkAnova"))
setGeneric("spkAccSD", function(object, spkSlopeOut, tol=3) standardGeneric("spkAccSD"))
setGeneric("spkPot", function(object, spkSlopeOut, sig, SD, precisionQuantile)
           standardGeneric("spkPot"))
