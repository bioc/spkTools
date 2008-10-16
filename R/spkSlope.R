## signal detect slope for the entire datasets
## and for each bin -- low, med, high
spkSlope <- function(object, label=NULL, cuts=c(.6,.99), ...){
              ## set palette
              palette(RColorBrewer:::brewer.pal(8,"Dark2"))
              tmp <- spkSplit(object)
              s <- tmp$s
              ns <- tmp$ns
              n <- spikeIn(s)
              e <- exprs(s)
              ens <- rowMeans(exprs(ns),na.rm=TRUE)
              ## regression of all the data
              Reg <- lm(as.vector(e)~as.vector(n))
              plot(x=jitter(as.vector(n)), y=as.vector(e),
                   xlab="Nominal Concentration",
                   ylab="Observed Expression Measure",
                   main=label, ...)
              ## compute means at each nom conc
              amt <- unique(as.vector(n))
              amt <- sort(amt[!is.na(amt)])
              avgE <- vector(length=length(amt))
              for(i in 1:length(amt)){
                  avgE[i] <- mean(e[n==amt[i]], na.rm=TRUE)
              }
            ## one bin
            if(is.null(cuts)){
              rSlps <- round(Reg$coef[2], digits=2)
              qout <- spkQuantile(amt, avgE, ens, p=c(0,1))
              brkpts <- c(-Inf,Inf)
              brks <- qout$out
              prop <- qout$prop
              ## plot slopes
              a <- Reg$coef[1]
              b <- Reg$coef[2]
              x0 <- min(amt[brks==min(brks)])
              x1 <- max(amt[brks==max(brks)])
              y0 <- a+b*x0
              y1 <- a+b*x1
              segments(x0, y0, x1, y1, col=3, lwd=4)
              ## legends
              txt <- paste("Slope:", rSlps[1])
              legend("bottomright", txt, col=3, lty=1, lwd=4, cex=.8,
                     bg="white")
              return(list(avgExp=avgE, slopes=rSlps, breaks=rbind(amt, brks),
                          brkpts=brkpts, prop=prop))
            }
            ## three bins  
            if(!is.null(cuts)){
              ## compute 3 slopes (low,med,high)
              qout <- spkQuantile(amt, avgE, ens, p=cuts)
              brkpts <- quantile(ens,probs=cuts,na.rm=TRUE)
              brks <- qout$out
              prop <- qout$prop
              B <- vector(length=3)
              lowInd <- n %in% amt[brks == 1]
              if(sum(brks == 1) > 1){
                  B[1] <- TRUE
                  lowE <- e[lowInd]
                  lowN <- n[lowInd]
                  lowReg <- lm(lowE~lowN)
              }
              medInd <- n %in% amt[brks == 2]
              if(sum(brks == 2) > 1){
                  B[2] <- TRUE
                  medE <- e[medInd]
                  medN <- n[medInd]
                  medReg <- lm(medE~medN)
              }
              highInd <- n %in% amt[brks == 3]
              if(sum(brks == 3) > 1){
                  B[3] <- TRUE
                  highE <- e[highInd]
                  highN <- n[highInd]
                  highReg <- lm(highE~highN)
              }
              ## plot slopes
              a <- Reg$coef[1]
              b <- Reg$coef[2]
              x0 <- min(amt[brks==min(brks)])
              x1 <- max(amt[brks==max(brks)])
              y0 <- a+b*x0
              y1 <- a+b*x1
              segments(x0, y0, x1, y1, col=3, lwd=4)
              if(B[1]){
                  a <- lowReg$coef[1]
                  b <- lowReg$coef[2]
                  x0 <- min(amt[brks==1])
                  x1 <- max(amt[brks==1])
                  y0 <- a+b*x0
                  y1 <- a+b*x1
                  segments(x0, y0, x1, y1, col=4, lwd=4)
              }
              if(B[2]){
                  a <- medReg$coef[1]
                  b <- medReg$coef[2]
                  x2 <- min(amt[brks==2])
                  x3 <- max(amt[brks==2])
                  y2 <- a+b*x2
                  y3 <- a+b*x3
                  segments(x2, y2, x3, y3, col=5, lwd=4)
              }
              if(B[3]){
                  a <- highReg$coef[1]
                  b <- highReg$coef[2]
                  x4 <- min(amt[brks==3])
                  x5 <- max(amt[brks==3])
                  y4 <- a+b*x4
                  y5 <- a+b*x5
                  segments(x4, y4, x5, y5, col=6, lwd=4)
              }

              ## add cut point lines
              if(B[1] & B[2]) abline(v=(x1+x2)/2, lty=2, lwd=3)
              if(B[2] & B[3]) abline(v=(x3+x4)/2, lty=2, lwd=3)

              rSlps <- rep(NA,4)
              rSlps[1] <- Reg$coef[2]
              if(B[1]) rSlps[2] <- lowReg$coef[2]
              if(B[2]) rSlps[3] <- medReg$coef[2]
              if(B[3]) rSlps[4] <- highReg$coef[2]
              rSlps <- round(rSlps, digits=2)
              names(rSlps) <- c("Total", "Low", "Med", "High")
              txt <- c(paste("Total:", rSlps[1]),
                       paste("Low:", rSlps[2]), paste("Med:", rSlps[3]),
                       paste("High:", rSlps[4]))
              legend("bottomright", txt, col=3:6, lty=1, lwd=4, cex=.8,
                     bg="white")
              return(list(avgExp=avgE, slopes=rSlps, breaks=rbind(amt, brks),
                          brkpts=brkpts, prop=prop))
            }
          }
          
