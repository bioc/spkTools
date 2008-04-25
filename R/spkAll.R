## run spike-in funcs (possibly for each sample)
setMethod("spkAll", "SpikeInExpressionSet",
          function(object, label, model=expr~spike+probe+array, fc=NULL, tol=3, compare=NULL, xrngs=NULL, yrngs=NULL, grps=NULL, cuts=c(.6,.99), potQuantile=.995, xaxis=NULL, pch=".", output="eps"){
              if("sample" %in% varLabels(object)){
                  ind <- pData(object)$sample
                  N <- max(ind)
                  for(i in 1:N){
                      tmp <- object[,ind == i]
                      spkFuncs(object=tmp, label=paste(arrayType, i, sep="."),
                               model, fc, tol, compare, xrngs, yrngs, grps,
                               cuts, potQuantile, xaxis, pch, output)
                  }
              } else spkFuncs(object, label, model, fc, tol, compare,
                              xrngs, yrngs, grps, cuts, potQuantile,
                              xaxis, pch, output)
          })





