## wrapper for spikein functions
spkAll <- function(object, label, model=expr~spike+probe+array, fc=NULL,
                   tol=3, xrngs=NULL, yrngs=NULL, cuts=c(.6,.99),
                   potQuantile=.995, gnn=c(25,100,10000), pch=".", output="eps"){
            mypar <- function(a=1,b=1,brewer.n=8,brewer.name="Dark2",...){
              par(mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
              par(mfrow=c(a,b),...)
              palette(RColorBrewer:::brewer.pal(brewer.n,brewer.name))
            }
            ## output label
            olab <- paste(strsplit(label, " ")[[1]], collapse="")
            
              ## Figure 1 - Slope Plot
              cat("Generating slope plot... ")
              if(output=="eps") postscript(paste("slope", olab, ".eps",
                   sep=""), horizontal=FALSE, encoding="TeXtext.enc",
                   pointsize=15, width=8, height=8)
              if(output=="pdf") pdf(paste("slope", olab, ".pdf", sep=""),
                   pointsize=15, width=8, height=8)
              mypar(mar=c(2.5,2.5,2,0.5), cex=1.8)
              ss <- spkSlope(object, label, cuts=cuts, pch=pch,
                             xlim=xrngs$s, ylim=yrngs$s)
              dev.off()
              cat("Done\n")

              ## Figure 2 - Density Plot
              cat("Generating density plot... ")
              if(output=="eps") postscript(paste("density", olab, ".eps",
                   sep=""), horizontal=FALSE, encoding="TeXtext.enc",
                   pointsize=15, width=8, height=8)
              if(output=="pdf") pdf(paste("density", olab, ".pdf", sep=""),
                   pointsize=15, width=8, height=8)
              mypar(mar=c(2.5,2.5,2,0.5), cex=1.8)
              drawcuts <- !is.null(cuts)
              spkDensity(object, ss, drawcuts, label,
                         xlim=xrngs$d, ylim=yrngs$d)
              dev.off()
              cat("Done\n")

              ## Table 1
              cat("Generating table 1... ")
              vtmp <- spkVar(object)
              sv <- as.numeric(vtmp[,2][-nrow(vtmp)])
              bin <- c("Low", "Med", "High")
              bins <- bin[ss$breaks[2,]]
              tab1 <- data.frame(NominalConc=2^ss$breaks[1,],
                                 AvgExp=round(ss$avgExp,1),
                                 PropGenesBelow=round(ss$prop,2),
                                 ALEStrata=bins,
                                 SD=round(sv,2))
              write.csv(tab1, file=paste("spkT1", olab, ".csv", sep=""),
                        row.names=FALSE)
              cat("Done\n")

              ## Figure 3 - Boxplots
              if(!is.null(fc) & !is.null(cuts)){
                  cat("Generating fold change plots... ")
                  spkBoxOut <- spkBox(object, spkSlopeOut=ss, fc=fc,
                                 tol=tol)
                  if(output=="eps") postscript(paste("FC", fc, "box", olab,
                       ".eps", sep=""), horizontal=FALSE,
                       encoding="TeXtext.enc", pointsize=15, width=8, height=8)
                  if(output=="pdf") pdf(paste("FC", fc, "box", olab, ".pdf",
                       sep=""), pointsize=15, width=8, height=8)
                  mypar(mar=c(5.2,2.5,2,0.5), cex=1.8)
                  plotSpkBox(spkBoxOut, fc=fc, xlim=xrngs$v, ylim=yrngs$v, main=label)
                  dev.off()
                  sbox <- summarySpkBox(spkBoxOut)
                  cat("Done\n")
                  
              ## Table 2
                  cat("Generating table 2... ")
                  AccuracySlope <- round(ss$slopes[-1], digits=2)
                  AccuracySD <- round(spkAccSD(object, ss), digits=2)
                  pot <- spkPot(object, ss, AccuracySlope, AccuracySD,
                                potQuantile)
                  PrecisionSD <- round(sbox$madFC[1:3], digits=2)
                  PrecisionQuantile <- round(pot$quantiles, digits=2)
                  SNR <- round(AccuracySlope/PrecisionSD, digits=2)
                  POT <- round(pot$POTs, digits=2)
                  GNN <- vector(length=3)
                  for(k in 1:3){
                    GNN[k] <- spkGNN(n=gnn[1], n.expr=gnn[2], n.unexpr=gnn[3],
                                     AccuracySlope[k], AccuracySD[k], spkBoxOut[[k]])
                  }
                  tab2 <- data.frame(AccuracySlope=AccuracySlope,
                                     AccuracySD=AccuracySD,
                                     PrecisionSD=PrecisionSD,
                                     PrecisionQuantile=PrecisionQuantile,
                                     SNR=SNR,
                                     POT=POT,
                                     GNN=GNN)
                  write.csv(tab2, file=paste("spkT2", olab, ".csv", sep=""),
                            row.names=c("Low","Med","High"))
                  cat("Done\n")
                  
              ## Figure 4 - MA plots
                  cat("Generating MA plots... ")
                  if(output=="eps") postscript(paste("FC", fc, "MA", olab,
                       ".eps", sep=""), horizontal=FALSE,
                       encoding="TeXtext.enc", pointsize=15, width=8, height=8)
                  if(output=="pdf") pdf(paste("FC", fc, "MA", olab, ".pdf",
                       sep=""), pointsize=15, width=8, height=8)
                  mypar(mar=c(5.2,2.5,2,0.5), cex=1.8)
                  spkMA(object, label=label, spkSlopeOut=ss, fc=fc,
                        tol=tol, ylim=yrngs$m)
                  dev.off()
                  cat("Done\n")
                }
            
              ## Limited Figure 3 & Table 2
              if(!is.null(fc) & is.null(cuts)){
                  cat("Generating fold change plots... ")
                  spkBoxOut <- spkBox(object, spkSlopeOut=ss, fc=fc,
                                 tol=tol)
                  if(output=="eps") postscript(paste("FC", fc, "box", olab,
                       ".eps", sep=""), horizontal=FALSE,
                       encoding="TeXtext.enc", pointsize=15, width=8, height=8)
                  if(output=="pdf") pdf(paste("FC", fc, "box", olab, ".pdf",
                       sep=""), pointsize=15, width=8, height=8)
                  mypar(mar=c(5.2,2.5,2,0.5), cex=1.8)
                  plotSpkBox(spkBoxOut, fc=fc, xlim=xrngs$v,
                              ylim=yrngs$v, main=label)
                  dev.off()
                  sbox <- summarySpkBox(spkBoxOut)
                  cat("Done\n")
                  
              ## Table 2
                  cat("Generating table 2... ")
                  AccuracySlope <- c(NA, round(ss$slopes, digits=2), NA)
                  AccuracySD <- round(spkAccSD(object, ss), digits=2)
                  pot <- spkPot(object, ss, AccuracySlope, AccuracySD,
                                potQuantile)
                  PrecisionSD <- round(sbox$madFC[1], digits=2)
                  PrecisionQuantile <- round(pot$quantiles, digits=2)
                  SNR <- round(AccuracySlope/PrecisionSD, digits=2)
                  POT <- round(pot$POTs, digits=2)
                  GNN <- spkGNN(n=gnn[1], n.expr=gnn[2], n.unexpr=gnn[3],
                                AccuracySlope[2], AccuracySD[2],
                                spkBoxOut[[2]])
                  tab2 <- data.frame(AccuracySlope=AccuracySlope[2],
                                     AccuracySD=AccuracySD[2],
                                     PrecisionSD=PrecisionSD,
                                     PrecisionQuantile=PrecisionQuantile[2],
                                     SNR=SNR[2],
                                     POT=POT[2],
                                     GNN=GNN)
                  write.csv(tab2, file=paste("spkT2", olab, ".csv", sep=""),
                            row.names=c("Total"))
                  cat("Done\n")

              ## Figure 4 - MA plots
                  cat("Generating MA plots... ")
                  if(output=="eps") postscript(paste("FC", fc, "MA", olab,
                       ".eps", sep=""), horizontal=FALSE,
                       encoding="TeXtext.enc", pointsize=15, width=8, height=8)
                  if(output=="pdf") pdf(paste("FC", fc, "MA", olab, ".pdf",
                       sep=""), pointsize=15, width=8, height=8)
                  mypar(mar=c(5.2,2.5,2,0.5), cex=1.8)
                  spkMA(object, label=label, spkSlopeOut=ss, fc=fc,
                        tol=tol, ylim=yrngs$m, reduce=FALSE)
                  dev.off()
                  cat("Done\n")
                }
            
              ## Table 3
              cat("Generating table 3... ")
              bals <- round(spkBal(object))
              anv <- round(spkAnova(object, model), digits=2)
              tab3 <- t(c(anv,bals))
              write.csv(tab3, file=paste("spkT3", olab, ".csv", sep=""),
                        row.names=FALSE)
              cat("Done\n")
          }






