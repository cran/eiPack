bayes.regress <- function(formula, data, sample = 1000, weights = NULL){
  fml <- as.formula(formula)
  D <- model.frame(fml, data = data)
  Y <- as.matrix(model.response(D))
  X <- model.matrix(formula, data = D)
  if (!is.null(weights)) {
    weights <- weights/sum(weights)
    w <- weights^(-1/2)
    X <- w * X
    Y <- w * Y
  }
  n <- nrow(X)
  nu <- n - ncol(X)
  s2 <- as.numeric((t(Y) %*% Y - t(Y) %*% X %*% solve(t(X) %*% X) %*% t(X) %*% Y)) / nu
  phi <- rchisq(sample, df = nu)
  sigma2 <- (nu * s2) / phi
  beta.hat <- solve(t(X) %*% X) %*% t(X) %*% Y
  tXX <- function(sigma2, X) sigma2 * solve(t(X) %*% X)
  Sigma <- lapply(sigma2, tXX, X)
  beta <- t(sapply(Sigma, mvrnorm, n = 1, mu = beta.hat))
  colnames(beta) <- colnames(X)
  beta
}
