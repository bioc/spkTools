spkPot <- function(object, spkSlopeOut, sig, SD, precisionQuantile=.995){
    nsM <- spkPairNS(object)
    gc()
    ## determine which background probes are L,M,H
    brkpts <- spkSlopeOut$brkpts
    eMeans <- rowMeans(exprs(spkSplit(object)$ns))
    gc()
    ## compute quantiles
    p9 <- vector(length=3)
    bgHigh <- as.vector(nsM[eMeans>=brkpts[2], ])
    p9[3] <- quantile(bgHigh,probs=precisionQuantile,na.rm=TRUE)
    rm(bgHigh)
    gc()
    bgMed <- as.vector(nsM[eMeans>=brkpts[1]&eMeans<brkpts[2], ])
    p9[2] <- quantile(bgMed,probs=precisionQuantile,na.rm=TRUE)
    rm(bgMed)
    gc()
    bgLow <- as.vector(nsM[eMeans<brkpts[1], ])
    rm(nsM)
    gc()
    p9[1] <- quantile(bgLow,probs=precisionQuantile,na.rm=TRUE)
    rm(bgLow)
    gc()
    ## compute POTs
    POTs <- vector(length=3)
    POTs[1] <- 1-pnorm(p9[1], mean=sig[1], sd=SD[1])
    POTs[2] <- 1-pnorm(p9[2], mean=sig[2], sd=SD[2])
    POTs[3] <- 1-pnorm(p9[3], mean=sig[3], sd=SD[3])
    return(data.frame(quantiles=p9, POTs=POTs))
}
