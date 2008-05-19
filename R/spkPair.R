## compute all pairwise array comparisons for spikein probes
spkPair <- function(object){
            tmp <- spkSplit(object)
            s <- tmp$s
            e <- exprs(s)
            n <- spikeIn(s)
            p <- permutations(n=ncol(n),r=2)
            mafc <- array(dim=c(nrow(n),nrow(p),4),dimnames=c("probes","arraypairs","M,A,N1,N2"))
            for(k in 1:nrow(p)){
              e1 <- e[,p[k,1]]
              e2 <- e[,p[k,2]]
              mafc[,k,1] <- e1-e2
              mafc[,k,2] <- (e1+e2)/2
              mafc[,k,3] <- n[,p[k,1]]
              mafc[,k,4] <- n[,p[k,2]]
            }
            rownames(mafc) <- rownames(e)
            return(mafc)
          }

