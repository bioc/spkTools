## function to determine empirical quantiles for the spikeins
spkQuantile <- function(amt, avgE, ens, p){
    obs <- as.vector(ens)
    N <- length(obs)
    M <- length(amt)
    q <- vector(length=M)
    out <- vector(length=2)
    for(i in 1:M){
        q[i] <- sum(obs<=avgE[i], na.rm=TRUE)/N
    }
    out[q <= p[1]] <- 1
    out[q <= p[2] & q > p[1]] <- 2
    out[q > p[2]] <- 3
    return(list(out=out, prop=q))
}
