
summary.I1 <- function(obj, ...) {
  p0    <- ncol(obj$Z0)
  p1    <- ncol(obj$Z1)
  p2    <- ifelse(is.null(obj$Z2), 0, ncol(obj$Z2))
  Time  <- nrow(obj$Z0)
  rank  <- ncol(obj$beta)
  S21   <- obj$M[-(1:(p0 + p1)), (p0 + 1):(p0 + p1), drop = F]
  S11   <- obj$M[(p0 + 1):(p0 + p1), (p0 + 1):(p0 + p1)]
  S12   <- obj$M[(p0 + 1):(p0 + p1), -(1:(p0 + p1))]
  S22   <- obj$M[-(1:(p0 + p1)), -(1:(p0 + p1)), drop = F]
  PI    <- if (ncol(obj$beta) != 0) tcrossprod(obj$alpha, obj$beta)
  OMega <- crossprod(obj$residuals) / Time
  
  #     SIgma <- if (!is.null(obj$Z2)) {
  #                 if (!is.null(PI)) { 
  #                   S22 - S21 %*% obj$beta %*% ginv(crossprod(obj$beta,S11) %*% obj$beta) %*% crossprod(obj$beta, S12)
  #                 }
  #                 else NULL
  #              }
  #              else NULL
  #     
  #     if (!is.null(obj$Z2) & !is.null(PI)) {
  #       se.PSi <- sqrt(diag(OMega) %o% diag(solve(SIgma)) / Time)
  #       t.PSi  <- obj$Psi / se.PSi
  #     }
  #     else se.PSi <- t.PSi <- NULL
  
  if (!is.null(obj$LongRunRes)) {
    hessList <- hessian(obj)
    se.aLpha <- hessList$seAlpha
    t.aLpha  <- obj$alpha / se.aLpha
    se.PI    <- hessList$sePi
    t.PI     <- (obj$alpha %*% t(obj$beta)) / se.PI
  } 
  else if (!is.null(PI)) {
    # Standard errors and t-statistics for alpha and PI without imposed restrictions
    se.aLpha <- sqrt(diag(OMega) %o% diag(ginv(crossprod(obj$beta,S11) %*% obj$beta)) / Time)
    t.aLpha	 <- obj$alpha / se.aLpha
    se.PI    <- sqrt(diag(OMega) %o% diag(obj$beta %*% tcrossprod(ginv(crossprod(obj$beta, S11) %*% obj$beta), obj$beta)) / Time)
    t.PI     <- PI / se.PI
  }
  else se.aLpha <- t.aLpha <- se.PI <- t.PI <- NULL
  
  # Get loglikelihood values
  ll <- -nrow(obj$Z0) * 0.5 * c(determinant(OMega)$modulus) - Time * p0 * 0.5 * (log(2 * pi) + 1)
  
  colnames(PI)
  tmp             <- obj
  tmp$beta        <- obj$beta
  tmp$alpha       <- obj$alpha
  tmp$t.alpha     <- t.aLpha
  tmp$Pi          <- PI
  tmp$t.Pi        <- t.PI
  # tmp$Psi         <- obj$Psi
  # tmp$t.Psi       <- t.PSi
  # tmp$Sigma       <- SIgma
  tmp$Omega       <- OMega
  tmp$eigenvalues <- obj$values
  tmp$logLik      <- ll
  tmp$rank        <- rank
  tmp$p           <- c(p0 = p0, p1 = p1, p2 = p2)
  tmp$lags        <- obj$lags
  tmp$residuals   <- obj$residuals
  
  class(tmp) <- 'summary.I1'
  
  return(tmp)
}

summary.I2 <- function(obj, ...) {
  p0 <- ncol(obj$Z0)
  p1 <- ncol(obj$Z1)
  p2 <- ncol(obj$Z2)
  p3 <- ifelse(is.null(obj$Z3), 0, ncol(obj$Z3))
  Time <- nrow(obj$Z0)
  obj$Omega <- crossprod(obj$residuals) / Time
  obj$logLik <- -nrow(obj$Z0) * 0.5 * c(determinant(obj$Omega)$modulus) -
    Time * p0 * 0.5 * log(2 * pi * exp(1))
  
  class(obj) <- 'summary.I2'
  
  return(obj)
}
