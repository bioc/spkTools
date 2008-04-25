## define our new class: SpikeInExpressionSet
require(Biobase)
setClass("SpikeInExpressionSet", contains="ExpressionSet")
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
