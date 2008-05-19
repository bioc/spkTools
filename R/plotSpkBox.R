## plot the results of the spkBox function
plotSpkBox <- function(boxs, fc=2, box.names=NULL, ...){
              par(las=2, mar=c(5.5,4,4,2))
              if(!is.null(box.names)) names(boxs) <- box.names
              boxplot(boxs, pch="-",
                      ylab="Observed Log-Ratio", ...)
              abline(h=0,lty=3)
              abline(h=log2(fc),lty=3)
}
