## plot a density curve of the non-spike-in genes and a multicolor rug for the
## density of the spike-in genes
setMethod("spkDensity", "SpikeInExpressionSet",
          function(object, spkSlopeOut, cuts=TRUE, label=NULL, ...){
              avgExp <- spkSlopeOut$avgExp
              tmp <- spkSplit(object)
              s <- tmp$s
              ns <- tmp$ns
              es <- exprs(s)
              nes <- exprs(ns)
              nsden <- density(nes, adjust=3, na.rm=TRUE)
              if(!"xlim"%in%objects()) xlim <- c(min(c(es, nes), na.rm=TRUE), max(c(es, nes), na.rm=TRUE))
              plot(nsden, xlab="Observed Expression Measure",
                   ylab="Density", main=label, ...)
              rug(avgExp)
           if(!is.null(cuts)){
              brkpts <- spkSlopeOut$brkpts
              abline(v=brkpts[1], lty=3)
              abline(v=brkpts[2], lty=3)
            }
          }
          )
