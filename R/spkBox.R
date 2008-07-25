## boxplots of each comparison at a given fc
## generate breaks by spkSlope(object)$breaks
spkBox <- function(object, spkSlopeOut, fc=2, tol=3, reduce=TRUE){
              nsM <- spkPairNS(object)
              gc()
              ## determine which background probes are L,M,H
              brkpts <- spkSlopeOut$brkpts
              eMeans <- rowMeans(exprs(spkSplit(object)$ns), na.rm=TRUE)
              if(brkpts[1] > -Inf){
                bgLow <- as.vector(nsM[eMeans<brkpts[1], ])
                if(reduce){
                  N <- length(bgLow)
                  M <- 5000
                  bgLow <- sort(bgLow)[seq(N/(2*M),N-N/(2*M),len=M)]
                }
                gc()
              } else bgLow <- numeric(0)
              bgMed <- as.vector(nsM[eMeans>=brkpts[1]&eMeans<brkpts[2], ])
              if(reduce){
                N <- length(bgMed)
                M <- 5000
                bgMed <- sort(bgMed)[seq(N/(2*M),N-N/(2*M),len=M)]
              }
              gc()
              if(brkpts[2] < Inf){
                bgHigh <- as.vector(nsM[eMeans>=brkpts[2], ])
                if(reduce){
                  N <- length(bgHigh)
                  M <- 5000
                  bgHigh <- sort(bgHigh)[seq(N/(2*M),N-N/(2*M),len=M)]
                }
                gc()
              } else bgHigh <- numeric(0)
              rm(nsM)
              gc()
              ## now for spike-ins
              mafc <- spkPair(object)
              lfc <- round(log2(fc), digits=tol)
              ind <- round(mafc[,,3]-mafc[,,4],digits=tol)==lfc
              ind2 <- round(mafc[,,3]-mafc[,,4],digits=tol)==0
              ind[is.na(ind)] <- ind2[is.na(ind2)] <- FALSE
              sM <- mafc[,,1][ind]
              nulM <- mafc[,,1][ind2]
              N1 <- mafc[,,3][ind]
              nulN1 <- mafc[,,3][ind2]
              N2 <- mafc[,,4][ind]
              nulN2 <- mafc[,,4][ind2]
              ## which nominal concs are in which bins
              breaks <- spkSlopeOut$breaks
              low <- breaks[1,][breaks[2,]==1]
              med <- breaks[1,][breaks[2,]==2]
              high <- breaks[1,][breaks[2,]==3]
              ## populate the boxs
              boxs <- list()
              boxs[[1]] <- c(bgLow, -bgLow)
              boxs[[2]] <- c(bgMed, -bgMed)
              boxs[[3]] <- c(bgHigh, -bgHigh)
              boxs[[4]] <- nulM[nulN1%in%low & nulN2%in%low]
              boxs[[5]] <- nulM[nulN1%in%med & nulN2%in%med]
              boxs[[6]] <- nulM[nulN1%in%high & nulN2%in%high]
              boxs[[7]] <- sM[N1%in%low & N2%in%low]
              boxs[[8]] <- sM[N1%in%med & N2%in%low]
              boxs[[9]] <- sM[N1%in%med & N2%in%med]
              boxs[[10]] <- sM[N1%in%high & N2%in%low]
              boxs[[11]] <- sM[N1%in%high & N2%in%med]
              boxs[[12]] <- sM[N1%in%high & N2%in%high]
              ## which bins to include
              ind <- lapply(boxs,length)==0
              posxnames <- c("Bg-Null LL", "Bg-Null MM", "Bg-Null HH",
                             "S-Null LL", "S-Null MM", "S-Null HH", "LL", "ML",
                             "MM", "HL", "HM", "HH")
              boxs[ind] <- NULL
              names(boxs) <- posxnames[!ind]
              return(boxs)
          }

