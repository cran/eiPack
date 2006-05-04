ei.MD.bayes <- function(formula, covariate = NULL, total = NULL, data,
                        lambda1 = 4, lambda2 = 2,
                        tune.list = NULL,
                        start.list = NULL,
                        sample = 1000,
                        thin = 1, burnin = 1000, verbose = 0,
                        ret.beta = 'r', ret.mcmc = TRUE,
                        usrfun = NULL){

if(!is.null(total)){
  if(!is.numeric(total)){
    if(is.character(total)){ total <- data[[total]]}else{
      total <- data[[deparse(substitute(total))]]}}}
  
if(is.null(covariate)){
  
  if(is.null(tune.list)){
    tune.alpha <- NULL
    tune.beta <- NULL}else{
      tune.alpha <- tune.list[[1]]
      tune.beta <- tune.list[[2]]}

  
  if(is.null(start.list)){
    start.alphas <- NULL
    start.betas <- NULL}else{
      start.alphas <- tune.list[[1]]
      start.betas <- tune.list[[2]]}
  
  if(is.null(usrfun)){output <- BayesMDei(formula, data, total=total, lambda1 = lambda1, lambda2 = lambda2,
                      tune.alpha = tune.alpha, tune.beta = tune.beta,
                                          start.alphas = start.alphas,
                                          start.betas = start.betas, sample = sample,
                      thin = thin, burnin = burnin, verbose = verbose, ret.beta =
                      ret.beta, ret.mcmc = ret.mcmc)
		output <- list(list(output$Alpha, output$Beta, output$cell.count), list(output$alpha.acc, 
output$beta.acc))
	names(output) <- c("draws", "acc.ratios")
	names(output$draws) <- c("Alpha", "Beta", "Cell.counts")
	names(output$acc.ratios) <- c("alpha.acc", "beta.acc")}
  else{output <- BayesMDei2(formula, data, total=total, lambda1 = lambda1, lambda2 = lambda2,
                      tune.alpha = tune.alpha, tune.beta = tune.beta,
                            start.alpha = start.alphas, start.betas =
                            start.betas,
                            sample = sample,
                      thin = thin, burnin = burnin, verbose = verbose, ret.beta =
                      ret.beta, ret.mcmc = ret.mcmc, usrfun = usrfun)
output <- list(list(output$Alpha, output$Beta, output$cell.count), list(output$alpha.acc,
output$beta.acc), output$usrfun) 
        names(output) <- c("draws", "acc.ratios","usrfun") 
        names(output$draws) <- c("Alpha", "Beta", "Cell.counts") 
        names(output$acc.ratios) <- c("alpha.acc", "beta.acc")
}
}else{

  if(is.null(tune.list)){
    tune.dr <- NULL
    tune.beta <- NULL
    tune.gamma <- NULL
    tune.delta <- NULL}else{
      tune.dr <- tune.list[[1]]
      tune.beta <- tune.list[[2]]
      tune.gamma <- tune.list[[3]]
      tune.delta <- tune.list[[4]]}


  
  if(is.null(start.list)){
    start.dr <- NULL
    start.betas <- NULL
    start.gamma <- NULL
    start.delta <- NULL}else{
      start.dr <- tune.list[[1]]
      start.betas <- tune.list[[2]]
      start.gamma <- tune.list[[3]]
      start.delta <- tune.list[[4]]}
  
  
  if(is.null(usrfun)){output <-  BayesMDei3cov(formula, covariate,
                                               total = total, data, lambda1 =
                                               lambda1, lambda2 = lambda2,
                      tune.dr = tune.dr, tune.beta = tune.beta,
                      tune.gamma=tune.gamma, tune.delta = tune.delta,
                                               start.dr = start.dr,
                                               start.betas =
                                               start.betas,
                                               start.gamma =
                                               start.gamma,
                                               start.delta =
                                               start.delta,
                                               sample = sample,
                      thin = thin, burnin = burnin, verbose = verbose, ret.beta =
                      ret.beta, ret.mcmc = ret.mcmc)

output <- list(list(output$Dr, output$Beta, output$Gamma, output$Delta, output$cell.count), list(output$dr.acc,   
output$beta.acc, output$gamma.acc))
        names(output) <- c("draws", "acc.ratios")
        names(output$draws) <- c("Dr", "Beta", "Gamma", "Delta","Cell.counts")
        names(output$acc.ratios) <- c("dr.acc", "beta.acc","gamma.acc")
}
    else{output <- BayesMDei4cov(formula, covariate, total = total, data, lambda1 =
                                 lambda1, lambda2 = lambda2,
                      tune.dr = tune.dr, tune.beta = tune.beta,
                      tune.gamma=tune.gamma, tune.delta = tune.delta,
                                 start.dr = start.dr,
                                 start.betas = start.betas,
                                 start.gamma = start.gamma,
                                 start.delta = start.delta,
                                 sample = sample,
                      thin = thin, burnin = burnin, verbose = verbose, ret.beta =
                      ret.beta, ret.mcmc = ret.mcmc, usrfun = usrfun)
output <- list(list(output$Dr, output$Beta, output$Gamma, output$Delta, output$cell.count), list(output$dr.acc,
output$beta.acc, output$gamma.acc), output$usrfun)
        names(output) <- c("draws", "acc.ratios","usrfun")
        names(output$draws) <- c("Dr", "Beta", "Gamma", "Delta","Cell.counts")
        names(output$acc.ratios) <- c("dr.acc", "beta.acc","gamma.acc")
}
}
  output$call <- match.call()
  class(output) <- "eiMD"
  return(output)
}
  



