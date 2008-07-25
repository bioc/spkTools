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
