## estimate the sd for spike-ins at the lowest possible fc in each bin
setMethod("spkAccSD","SpikeInExpressionSet",
          function(object, spkSlopeOut, tol=3){
              breaks <- spkSlopeOut$breaks
              n <- m <- vector(length=3)
              for(k in 1:3){
                  tmp <- breaks[1,breaks[2,]==k]
                  n[k] <- length(tmp)
                  if (n[k] == 0) m[k] <- NA else m[k] <- round(min(tmp[2:n[k]]-tmp[1:(n[k]-1)]), tol)
              }

              um <- unique(round(m, tol))
              um <- um[!is.na(um)]
              N <- length(um)
              MADs <- matrix(nrow=3, ncol=N)
              for(k in 1:N){
                  fc <- round(um[k], digits=tol)
                  mafc <- spkPair(object)
                  ind <- round(mafc[,,3]-mafc[,,4], digits=tol)==fc
                  ind[is.na(ind)] <- FALSE
                  sM <- mafc[,,1][ind]
                  N1 <- mafc[,,3][ind]
                  N2 <- mafc[,,4][ind]
                  ## which nominal concs are in which bins
                  low <- breaks[1,][breaks[2,]==1]
                  med <- breaks[1,][breaks[2,]==2]
                  high <- breaks[1,][breaks[2,]==3]
                  ## populate the boxs
                  boxs <- list()
                  boxs[[1]] <- sM[N1%in%low & N2%in%low]
                  boxs[[2]] <- sM[N1%in%med & N2%in%med]
                  boxs[[3]] <- sM[N1%in%high & N2%in%high]
                  ind <- lapply(boxs,length)==0
                  MADs[,k] <- sapply(boxs, function(x){if (!prod(is.na(x))) mad(x, na.rm=TRUE) else NA})
              }
              madFC <- vector(length=3)
              names(madFC) <- c("LL","MM","HH")
              for(k in 1:3){
                if (is.na(m[k])) madFC[k] <- NA else madFC[k] <- MADs[k,which(um==m[k])]
              }
              return(madFC)
              })
