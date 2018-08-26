summary.I1 <- function(object, ...) {
  p0 <- ncol(object$Z0)
  p1 <- ncol(object$Z1)
  p2 <- ifelse(is.null(object$Z2), 0, ncol(object$Z2))
  Time <- nrow(object$Z0)
  rank <- ncol(object$beta)
  
  PI <- if (ncol(object$beta) != 0) tcrossprod(object$alpha, object$beta)
  
  OMega <- crossprod(object$residuals) / Time
  SIgma <- if (!is.null(object$Z2)) {
    if (!is.null(PI)) object$M[-(1:(p0 + p1)), -(1:(p0 + p1)), drop = F] -
      object$M[-(1:(p0 + p1)), (p0 + 1):(p0 + p1), drop = F] %*%
      object$beta %*%
      ginv(crossprod(object$beta,
                     object$M[(p0 + 1):(p0 + p1), (p0 + 1):(p0 + p1)]) %*%
             object$beta
      ) %*%
      crossprod(object$beta,
                object$M[(p0 + 1):(p0 + p1), -(1:(p0 + p1))]
      )
  }
  else NULL
  
  if (!is.null(PI)) {
    se.aLpha <- sqrt(diag(OMega) %o%
                       diag(ginv(crossprod(object$beta,
                                           object$S11) %*% object$beta)) / Time
    )
    t.aLpha	<- object$alpha / se.aLpha
    
    se.PI <- sqrt(diag(OMega) %o%
                    diag(object$beta %*%
                           tcrossprod(ginv(crossprod(object$beta, object$S11) %*%
                                             object$beta), object$beta)) / Time
    )
    t.PI <- PI / se.PI
  }
  else se.aLpha <- t.aLpha <- se.PI <- t.PI <- NULL
  
  ll <- -nrow(object$Z0) * 0.5 * c(determinant(OMega)$modulus) -
    Time * p0 * 0.5 * (log(2 * pi) + 1)
  
  if (!is.null(object$Z2) & !is.null(PI)) {
    se.PSi <- sqrt(diag(OMega) %o% diag(solve(SIgma)) / Time)
    t.PSi <- object$Psi / se.PSi
  }
  else se.PSi <- t.PSi <- NULL
  colnames(PI)
  
  tmp	<- object['call']
  tmp$beta <- object$beta
  tmp$alpha <- object$alpha
  tmp$t.alpha <- t.aLpha
  tmp$Pi <- PI
  tmp$t.Pi <- t.PI
  tmp$Psi <- object$Psi
  tmp$t.Psi <- t.PSi
  tmp$Omega <- OMega
  tmp$Sigma <- SIgma
  tmp$eigenvalues <-object$values
  tmp$logLik <- ll
  tmp$rank <- rank
  tmp$p <- c(p0 = p0, p1 = p1, p2 = p2)
  tmp$lags <- object$lags
  tmp$residuals <- object$residuals
  
  class(tmp) <- 'summary.I1'
  
  return(tmp)
}

summary.I2 <- function(object, ...) {
  p0 <- ncol(object$Z0)
  p1 <- ncol(object$Z1)
  p2 <- ncol(object$Z2)
  p3 <- ifelse(is.null(object$Z3), 0, ncol(object$Z3))
  Time <- nrow(object$Z0)
  object$Omega <- crossprod(object$residuals) / Time
  object$logLik <- -nrow(object$Z0) * 0.5 * c(determinant(object$Omega)$modulus) -
    Time * p0 * 0.5 * log(2 * pi * exp(1))
  
  class(object) <- 'summary.I2'
  
  return(object)
}
