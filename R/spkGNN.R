spkGNN <- function(n, n.expr, n.unexpr, AccuracySlope, AccuracySD, nullfc){
  nullfc <- nullfc[!is.na(nullfc)]
  
  N <- vector(length=10000)
  for(i in 1:10000){
    tmp <- sample(nullfc, n.unexpr, replace=TRUE)
    x <- qnorm(1-(n/n.expr), mean=AccuracySlope, sd=AccuracySD)
    n0 <- sum(tmp>x)
    N[i] <- n + n0
  }

  return(ceiling(mean(N)))
}
