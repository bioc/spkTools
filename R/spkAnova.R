
## anova: expr ~ spike + array + probe + error
spkAnova <- function(object, model=expr~spike+probe+array){
              ## check model
              chk <- terms.formula(model)
              covs <- labels(chk)
              Y <- rownames(attr(chk, "factors"))[1]
              tst <- prod(covs %in% c("spike","probe","array"))
              tst2 <- (Y == "expr")
              if((tst*tst2) != TRUE) stop("Invalid model specification.")

              s <- spkSplit(object)$s
              n <- spikeIn(s)
              e <- exprs(s)
              ## computes number of replicate arrays
              ## NOTE: number of replicates must be a constant!
              ifelse("array"%in%covs, reps <- 1,
                      reps <- ncol(n)/sum(!duplicated(n,MARGIN=2)))

              ## data for lm fit
              nr <- nrow(e)
              nc <- ncol(e)
              dat <- data.frame(expr=as.vector(e), spike=as.factor(as.vector(n)),
                                array=as.factor(rep(1:nc, rep(nr,nc))), probe=as.factor(rep(1:nr,nc)))
              ## fit model
              options(contrasts=c("contr.sum","contr.sum"))
              fit <- lm(model, data=dat)

              ## compute SS
              ssr <- sum(fit$residuals^2)
              dfr <- df.residual(fit)
              comp <- fit$coef
              asgn <- fit$assign
              plus1 <- -unlist(lapply(split(comp, asgn), sum))[-1]
              comp <- c(comp,plus1)
              asgn <- c(asgn,1:length(plus1))
              ss <- c(unlist(lapply(split(comp^2, asgn), sum)), ssr)[-1]
              df <- c(unlist(lapply(split(asgn, asgn), length)), dfr)[-1]
              ms <- ss/(df+1)
              i <- length(ms)

              out <- sqrt(c(ms[-i]/reps,ms[i]))
              names(out) <- c(covs,"error")

              return(out)
          }

