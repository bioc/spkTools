## summarize the results of the spkBox function
summarySpkBox <- function(boxs){
              meanFC <- sapply(boxs, mean, na.rm=TRUE)
              madFC <- sapply(boxs, function(x){if (!prod(is.na(x))) mad(x, na.rm=TRUE) else NA})
              return(data.frame(meanFC, madFC))
}
