## quantification of the imbalance of the spike-in design
setMethod("spkBal","SpikeInExpressionSet",
          function(object){
              s <- spkSplit(object)$s
              n <- spikeIn(s)
              e <- exprs(s)
              n[is.na(n)] <- -99
              trts <- sort(unique(as.vector(n)))
              T <- length(trts)
              ## number of probes (U1) and arrays (U2)
              U1 <- nrow(n)
              U2 <- ncol(n)

              ## PROBE IMBALANCE:
              pExp <- U2^2 / T
              prbs <- apply(n,1,function(x) sum(table(x)^2))
              pBal <- sum(prbs-pExp) / U1

              ## ARRAY IMBALANCE:
              aExp <- U1^2 / T
              prbs <- apply(n,2,function(x) sum(table(x)^2))
              aBal <- sum(prbs-aExp) / U2

              out <- c(pBal,aBal)
              names(out) <- c("Probe Imbalance", "Array Imbalance")
              return(out)
          })
