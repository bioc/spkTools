## generic functions
setGeneric("spikeIn", function(object) standardGeneric("spikeIn"))
setGeneric("spkSplit", function(object) standardGeneric("spkSplit"))
setGeneric("spikeIn<-", function(object, value) standardGeneric("spikeIn<-"))
setGeneric("createSpikeInExpressionSet",
           function(object, exprs, spikeIn, ...) standardGeneric("createSpikeInExpressionSet"))
