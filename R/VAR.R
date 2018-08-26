VAR <- function(object, ...) UseMethod('VAR')

VAR.I1 <- function(object, ...) {
  p <- ncol(object$X)
  
  A <- array(0, c(p, p, object$lags))
  A[,,1] <- with(object, diag(p) + tcrossprod(alpha, beta[seq_len(p),, drop = FALSE]))
  for (i in seq_len(object$lags - 1)) {
    A[, , i] <- A[, , i] + object$Psi[, seq_len(p) + (i - 1) * p]
    A[, , i + 1] <- A[, , i + 1] - object$Psi[, seq_len(p) + (i - 1) * p]
  }
  
  class(A) <- 'VAR'
  return(A)
}

VAR.I2 <- function(object, ...) {
  p <- ncol(object$X)
  
  A <- array(0, c(p, p, object$lags))
  
  A[, , 1] <- with(object, 2 * diag(p) +
                     tcrossprod(alpha, tau[seq_len(p), ] %*% rho + delta[seq_len(p),]) +
                     tcrossprod(zeta, tau[seq_len(p), ]))
  A[, , 2] <- with(object, -(diag(p) + tcrossprod(alpha, delta[seq_len(p), ]) +
                               tcrossprod(zeta, tau[seq_len(p), ])))
  
  for (i in seq_len(max(object$lags - 2, 0))) {
    A[, , i] <- A[, , i] + object$Psi[, seq_len(p) + (i - 1) * p]
    A[, , i + 1] <- A[, , i + 1] - 2 * object$Psi[, seq_len(p) + (i - 1) * p]
    A[, , i + 2] <- A[, , i + 2] + object$Psi[, seq_len(p) + (i - 1) * p]
  }
  
  class(A) <- 'VAR'
  return(A)
}
