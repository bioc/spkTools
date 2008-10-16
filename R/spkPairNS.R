## compute all pairwise array comparisons for 1000 non-spikein probes
## only returns the positive fold change for each pairing
spkPairNS <- function(object,output="M"){
            tmp <- spkSplit(object)
            ns <- tmp$ns
            e <- exprs(ns)
            p <- gtools:::combinations(n=ncol(e),r=2)
            m <- matrix(nrow=nrow(e),ncol=nrow(p))
            ## rows are probes
            ## cols are arraypairs
            if(output=="M"){
              for(k in 1:nrow(p)){
                m[,k] <- abs(e[,p[k,1]]-e[,p[k,2]])
              }
            }
            if(output=="A"){
              for(k in 1:nrow(p)){
                m[,k] <- (e[,p[k,1]]+e[,p[k,2]])/2
              }
            } 
            return(m)
          }
          

