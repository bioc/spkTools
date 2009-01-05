setMethod("initialize", "SpikeInExpressionSet",
          function(.Object,
                   assayData = assayDataNew(exprs=exprs, spikeIn=spikeIn, ...),
                   exprs=new("matrix"),
                   spikeIn=new("matrix"),
                   phenoData = annotatedDataFrameFrom  (assayData, byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData = new("MIAME"),
                   annotation = new("character"), ...){
            .Object <- callNextMethod(.Object,
                                      assayData = assayData,
                                      phenoData = phenoData,
                                      experimentData = experimentData,
                                      annotation = annotation,
                                      featureData = featureData)
            .Object
          })

setMethod("createSpikeInExpressionSet",
          signature(object="SpikeInExpressionSet", exprs="matrix", spikeIn="matrix"),
          function(object, exprs, spikeIn, ...) new("SpikeInExpressionSet", exprs=exprs, spikeIn=spikeIn, ...))


setMethod("spikeIn", signature(object="SpikeInExpressionSet"),
          function(object) assayDataElement(object, "spikeIn")
          )

setReplaceMethod("spikeIn", signature(object="SpikeInExpressionSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "spikeIn", value))

## function to make a new SpikeInExpressionSet object containing only
## the spiked in genes and one containing only the non-spiked in genes
setMethod("spkSplit", "SpikeInExpressionSet",
          function(object){
              ## genes spikedin on any of the arrays
              spk <- as.logical(apply(!is.na(spikeIn(object)), 1, sum,
                                      na.rm=TRUE))
              s <- new("SpikeInExpressionSet", exprs=exprs(object)[spk,],
                       spikeIn=spikeIn(object)[spk,])
              ns <- new("SpikeInExpressionSet", exprs=exprs(object)[!spk,],
                        spikeIn=spikeIn(object)[!spk,])
              list(s=s, ns=ns)
          }
          )
