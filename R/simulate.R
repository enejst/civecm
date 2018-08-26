simulate.I1 <- function(object, nsim, seed, res, init, ...) {
  stopifnot(!missing(res))
  
  p <- nrow(object$alpha)
  r <- ncol(object$alpha)
  k <- object$lags
  
  A <- VAR(object)
  
  if(missing(init)) init <- object$X[seq_len(k), ]
  
  eps <- as.matrix(res) +
    tcrossprod(object$Z1[, seq_len(ncol(object$Z1) - p) + p, drop = F] %*%
                 object$beta[seq_len(ncol(object$Z1) - p) + p, , drop = F],
               object$alpha
    )
  if (!is.null(ncol(object$Z2)) && ncol(object$Z2) > (k - 1) * p)
    eps <-  eps + tcrossprod(object$Z2[, seq_len(ncol(object$Z2) - (k - 1) * p) +
                                         (k - 1) * p, drop = F],
                             object$Psi[, seq_len(ncol(object$Z2) - (k - 1) * p) +
                                          (k - 1) * p, drop = F]
    )
  attributes(eps) <- attributes(res)
  ans <- rbind(init, rmfilter(eps, A, init))
  
  return(ans)
}

simulate.I2 <- function(object, nsim, seed, res, init, ...) {
  stopifnot(!missing(res))
  
  p <- nrow(object$alpha)
  r <- ncol(object$alpha)
  k <- object$lags
  
  A <- VAR(object)
  
  if(missing(init)) init <- object$X[seq_len(k), ]
  
  eps <- as.matrix(res) +
    tcrossprod(object$Z2[, seq_len(ncol(object$Z2) - p) + p, drop = F] %*%
                 object$beta[seq_len(ncol(object$Z2) - p) + p, , drop = F] +
                 object$Z1[, seq_len(ncol(object$Z1) - p) + p, drop = F] %*%
                 object$delta[seq_len(ncol(object$Z1) - p) + p, , drop = F],
               object$alpha
    ) + object$Z1[, seq_len(ncol(object$Z1) - p) + p, drop = F] %*%
    tcrossprod(object$tau[seq_len(ncol(object$Z1) - p) + p,
                          , drop = F],
               object$zeta
    )
  if (!is.null(ncol(object$Z3)) && ncol(object$Z3) > (k - 2) * p)
    eps <- eps + tcrossprod(object$Z3[, seq_len(ncol(object$Z3) - (k - 2) * p) +
                                        (k - 2) * p, drop = F],
                            object$Psi[, seq_len(ncol(object$Z3) - (k - 2) * p) +
                                         (k - 2) * p, drop = F]
    )
  attributes(eps) <- attributes(res)
  
  ans <- rbind(init, rmfilter(eps, A, init))
  
  return(ans)
}
