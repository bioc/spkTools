## variance of expression at each spikein level
setMethod("spkVar", "SpikeInExpressionSet",
          function(object){
              tmp <- spkSplit(object)
              s <- tmp$s
              ns <- tmp$ns
              n <- spikeIn(s)
              e <- exprs(s)
              ens <- exprs(ns)
              lvls <- sort(unique(as.vector(n)))
              K <- length(lvls)
              V <- vector(length=K+1)
              for(k in 1:K){
                  ind <- n == lvls[k]
                  ind[is.na(ind)] <- FALSE
                  V[k] <- mad(as.vector(e[ind]), na.rm=TRUE)
              }
              V[K+1] <- mad(as.vector(ens), na.rm=TRUE)
              V <- cbind(c(lvls, "NonSpikeIns"), V)
              colnames(V) <- c("SpikeIn Amount", "MAD")
              return(V)
          }
          )
