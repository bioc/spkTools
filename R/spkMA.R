## boxplots of each comparison at a given fc
## generate breaks by spkSlope(object)$breaks
spkMA <- function(object, spkSlopeOut, fc=2, tol=3, label=NULL, ylim=NULL,
                   outlier=1, reduce=TRUE, plot.legend=TRUE){
              ## M: log ratio
              nsM <- spkPairNS(object,output="M")
              gc()
              ## determine which background probes are L,M,H
              brkpts <- spkSlopeOut$brkpts
              eMeans <- rowMeans(exprs(spkSplit(object)$ns), na.rm=TRUE)
              bgLowM <- as.vector(nsM[eMeans<brkpts[1], ])
              if(reduce){
                N <- length(bgLowM)
                M <- 5000
                bgLowMsort <- sort(bgLowM, index.return=TRUE)
                bgLowInd <- bgLowMsort$ix
                bgLowM <- bgLowMsort$x[seq(N/(2*M),N-N/(2*M),len=M)]
              }
              gc()
              bgMedM <- as.vector(nsM[eMeans>=brkpts[1]&eMeans<brkpts[2], ])
              if(reduce){
                N <- length(bgMedM)
                M <- 5000
                bgMedMsort <- sort(bgMedM, index.return=TRUE)
                bgMedInd <- bgMedMsort$ix
                bgMedM <- bgMedMsort$x[seq(N/(2*M),N-N/(2*M),len=M)]
              }
              gc()
              bgHighM <- as.vector(nsM[eMeans>=brkpts[2], ])
              if(reduce){
                N <- length(bgHighM)
                M <- 5000
                bgHighMsort <- sort(bgHighM, index.return=TRUE)
                bgHighInd <- bgHighMsort$ix
                bgHighM <- bgHighMsort$x[seq(N/(2*M),N-N/(2*M),len=M)]
              }
              gc()
              rm(nsM)
              gc()
              ## A: avg exprs
              nsA <- spkPairNS(object,output="A")
              gc()
              ## determine which background probes are L,M,H
              brkpts <- spkSlopeOut$brkpts
              eMeans <- rowMeans(exprs(spkSplit(object)$ns), na.rm=TRUE)
              bgLowA <- as.vector(nsA[eMeans<brkpts[1], ])
              if(reduce){
                N <- length(bgLowA)
                M <- 5000
                bgLowA <- bgLowA[bgLowInd][seq(N/(2*M),N-N/(2*M),len=M)]
              }
              gc()
              bgMedA <- as.vector(nsA[eMeans>=brkpts[1]&eMeans<brkpts[2], ])
              if(reduce){
                N <- length(bgMedA)
                M <- 5000
                bgMedA <- bgMedA[bgMedInd][seq(N/(2*M),N-N/(2*M),len=M)]
              }
              gc()
              bgHighA <- as.vector(nsA[eMeans>=brkpts[2], ])
              if(reduce){
                N <- length(bgHighA)
                M <- 5000
                bgHighA <- bgHighA[bgHighInd][seq(N/(2*M),N-N/(2*M),len=M)]
              }
              gc()
              rm(nsA)
              gc()
              ## now for spike-ins
              mafc <- spkPair(object)
              lfc <- round(log2(fc), digits=tol)
              ind <- round(mafc[,,3]-mafc[,,4],digits=tol)==lfc
              ind2 <- round(mafc[,,3]-mafc[,,4],digits=tol)==0
              ind[is.na(ind)] <- ind2[is.na(ind2)] <- FALSE
              sM <- mafc[,,1][ind]
              sA <- mafc[,,2][ind]
              nulM <- mafc[,,1][ind2]
              nulA <- mafc[,,2][ind2]
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
              M <- A <- list()
              M[[1]] <- c(bgLowM, -bgLowM)
              M[[2]] <- c(bgMedM, -bgMedM)
              M[[3]] <- c(bgHighM, -bgHighM)
              M[[4]] <- nulM[nulN1%in%low & nulN2%in%low]
              M[[5]] <- nulM[nulN1%in%med & nulN2%in%med]
              M[[6]] <- nulM[nulN1%in%high & nulN2%in%high]
              M[[7]] <- sM[N1%in%low & N2%in%low]
              M[[8]] <- sM[N1%in%med & N2%in%low]
              M[[9]] <- sM[N1%in%med & N2%in%med]
              M[[10]] <- sM[N1%in%high & N2%in%low]
              M[[11]] <- sM[N1%in%high & N2%in%med]
              M[[12]] <- sM[N1%in%high & N2%in%high]
              A[[1]] <- c(bgLowA, bgLowA)
              A[[2]] <- c(bgMedA, bgMedA)
              A[[3]] <- c(bgHighA, bgHighA)
              A[[4]] <- nulA[nulN1%in%low & nulN2%in%low]
              A[[5]] <- nulA[nulN1%in%med & nulN2%in%med]
              A[[6]] <- nulA[nulN1%in%high & nulN2%in%high]
              A[[7]] <- sA[N1%in%low & N2%in%low]
              A[[8]] <- sA[N1%in%med & N2%in%low]
              A[[9]] <- sA[N1%in%med & N2%in%med]
              A[[10]] <- sA[N1%in%high & N2%in%low]
              A[[11]] <- sA[N1%in%high & N2%in%med]
              A[[12]] <- sA[N1%in%high & N2%in%high]
              ## all background
              bgM <- c(M[[1]],M[[2]],M[[3]])
              bgA <- c(A[[1]],A[[2]],A[[3]])
              ind <- lapply(M,length)>0
              posxnames <- c("Bg-Null LL", "Bg-Null MM", "Bg-Null HH",
                             "S-Null LL", "S-Null MM", "S-Null HH",
                             "LL", "ML", "MM", "HL", "HM", "HH")
              xaxis <- posxnames[ind]
              i.spk <- which(xaxis%in%c("S-Null LL", "S-Null MM", "S-Null HH"))
              legend.names <- xaxis[-i.spk]
              nulpch <- 22
              spkpch <- 17
              spkcols <- brewer.pal(8,"Set1")[c(1,3,5:8)]
              nulcols <- brewer.pal(8,"Greys")[6:8]
              ifelse(is.null(label), ptitle <- paste("FC =", fc),
                     ptitle <- label)
              iOutlier <- bgM>=outlier | bgM<=-outlier
              if(is.null(ylim)){
                ylim <- c(min(unlist(M)),max(unlist(M)))
                smoothScatter(y=bgM[!iOutlier], x=bgA[!iOutlier],
                              main=ptitle, nrpoints=0,
                              xlab="A", ylab="M", ylim=ylim)
              } else{
                smoothScatter(y=bgM[!iOutlier], x=bgA[!iOutlier],
                              main=ptitle, nrpoints=0,
                              xlab="A", ylab="M", ylim=ylim)
              }
              points(y=bgM[iOutlier],x=bgA[iOutlier],pch=".",col="blue",cex=2)
              ## for(k in 1:3){
              ##   if(length(M[k])>0) points(y=M[[k+3]],x=A[[k+3]],pch=nulpch,col=nulcols[k],cex=.4)
              ## }
              for(k in 1:6){
                if(length(M[k])>0) points(y=M[[k+6]],x=A[[k+6]],pch=spkpch,col=spkcols[k],cex=.5)
              }
              ## pchs <- c(rep(nulpch,sum(ind[4:6])),rep(spkpch,sum(ind[7:12])))
              ## cols <- c(nulcols,spkcols)[ind[4:12]]
              pchs <- rep(spkpch,sum(ind[7:12]))
              cols <- spkcols[ind[7:12]]
              par(lwd=2,cex=1.5)
              ## determine how many bgs are present
              bgs <- vector(length=3)
              bgs[1] <- length(M[[1]])>0
              bgs[2] <- length(M[[2]])>0
              bgs[3] <- length(M[[3]])>0
              n.bg <- sum(bgs)              
              if(plot.legend){
                legend("bottomright", legend.names, pch=c(rep(15,n.bg),pchs),
                       col=c(rep("blue",n.bg),cols), ncol=2, cex=.8,
                       pt.cex=.8, box.lwd=1)
              }
          }
