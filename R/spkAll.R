## run spike-in funcs (possibly for each sample)
spkAll <- function(object, label, model=expr~spike+probe+array, fc=NULL, tol=3,
                   xrngs=NULL, yrngs=NULL, cuts=c(.6,.99), potQuantile=.995,
                   pch=".", output="eps"){
              if("sample" %in% varLabels(object)){
                  ind <- pData(object)$sample
                  N <- max(ind)
                  for(i in 1:N){
                      tmp <- object[,ind == i]
                      spkFuncs(object=tmp, label=paste(arrayType, i, sep="."),
                               model, fc, tol, xrngs, yrngs,
                               cuts, potQuantile, pch, output)
                  }
              } else spkFuncs(object, label, model, fc, tol, 
                              xrngs, yrngs, cuts, potQuantile,
                              pch, output)
          }





