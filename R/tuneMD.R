tuneMD <- function(formula, covariate = NULL, data, ntunes = 10, sample = 1000, 
...) {
  D <- model.frame(formula, data)
  Groups <- D[[2]]
  ng <- ncol(Groups)
  N <- t(D[[1]])
  np <- nrow(N)
  npm1 <- np-1
  precincts <- nrow(D)
  sdtune <- function(xx){
    if((.3 < xx) & (xx <= .45)) return(.95)
    if(xx <= .3) return(.84)
    if((.45 < xx) & (xx < .7)) return(1.05)
    if(.7 <= xx) return(1.15)
  }
  
  if (is.null(covariate)) { 
    tuneA <- matrix(0.25, nrow = ng, ncol = np)
    tuneB <- array(0.05, dim = c(ng, npm1, precincts))
    
    for(jj in 1:ntunes){
      tl <- list(tune.alpha = tuneA, tune.beta = tuneB)
      tmp <- ei.MD.bayes(formula, covariate = covariate, data = data, sample = sample,
                         tune = tl, ...)
      Beta <- array(tmp$acc.ratios$beta.acc, dim = c(ng, np, precincts))
      Alpha <- matrix(tmp$acc.ratios$alpha.acc, nrow = ng, ncol = np)
      for (ii in 1:precincts) {
        for(rr in 1:ng){
          for(cc in 1:npm1){
            tuneB[rr,cc,ii] <- tuneB[rr,cc,ii] * sdtune(Beta[rr,cc,ii])
          }
        }
      }
      for(rr in 1:ng){
        for(cc in 1:np){
          tuneA[rr,cc] <- tuneA[rr,cc] * sdtune(Alpha[rr,cc])
        }
      }
    }
    return(list(tune.alpha = tuneA, tune.beta = tuneB))
  }
  else {
    tuneDr <- array(0.05, ng)
    tuneB <- array(0.05, dim = c(ng, npm1, precincts))
    tuneD <- tuneG <-  matrix(0.25, nrow = ng, ncol = npm1)
    for(jj in 1:ntunes){
      tl <- list(tune.dr = tuneDr, tune.beta = tuneB, tune.gamma = tuneG,
                 tune.delta = tuneD)
      tmp <- ei.MD.bayes(formula, covariate = covariate, data = data, sample = sample,
                         tune.list = tl, ...)
      Dr <- tmp$acc.ratios$dr.acc
      Beta <- array(tmp$acc.ratios$beta.acc, dim = c(ng, np, precincts))
      Gamma <- matrix(tmp$acc.ratios$gamma.acc, nrow = ng, ncol = npm1)
      Delta <- matrix(tmp$acc.ratios$gamma.acc, nrow = ng, ncol = npm1)
      for (ii in 1:precincts) {
        for(rr in 1:ng){
          for(cc in 1:npm1){
            tuneB[rr,cc,ii] <- tuneB[rr,cc,ii] * sdtune(Beta[rr,cc,ii])
          }
        }
      }
      for(rr in 1:ng){
        tuneDr[rr] <- tuneDr[rr] * sdtune(Dr[rr])
        for(cc in 1:npm1){
          tuneG[rr,cc] <- tuneG[rr,cc]*sdtune(Gamma[rr,cc])
          tuneD[rr,cc] <- tuneD[rr,cc]*sdtune(Delta[rr,cc])
        }
      }
    }
    return(list(tune.dr = tuneDr, tune.beta = tuneB, tune.gamma = tuneG, 
                tune.delta = tuneD))
  } 
}
