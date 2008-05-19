## multicolor rug function
multiColRug <- function(den, col="YlOrRd", ticksize=0.03, side=1, lwd=0.5) 
{
    x <- den$x
    y <- den$y
    n <- length(y)
    ind <- rank(y)
    rampcol <- colorRampPalette(brewer.pal(9, col))(n)
    for(i in 1:n){
      Axis(side=side, at=x[i], labels=FALSE, lwd=lwd, col=rampcol[ind[i]],
           tck=ticksize)
    }
}
