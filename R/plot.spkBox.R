## plot the results of the spkBox function
plot.spkBox <- function(boxs, fc=2, ...){
              par(las=2, cex=1.1)
              boxplot(boxs, pch="-",
                      ylab="Observed Log-Ratio", ...)
              abline(h=0,lty=3)
              abline(h=log2(fc),lty=3)
}
