simulate.I1 <- function(obj, nsim, seed, res, init, ...) {
  stopifnot(!missing(res))
  
  p <- nrow(obj$alpha)
  r <- ncol(obj$alpha)
  k <- obj$lags
  
  A <- VAR(obj)
  
  if(missing(init)) init <- obj$X[seq_len(k), ]
  
  eps <- as.matrix(res) +
    tcrossprod(obj$Z1[, seq_len(ncol(obj$Z1) - p) + p, drop = F] %*%
                 obj$beta[seq_len(ncol(obj$Z1) - p) + p, , drop = F],
               obj$alpha
    )
  
  if (!is.null(ncol(obj$Z2)) && ncol(obj$Z2) > (k - 1) * p)
    eps <-  eps + tcrossprod(obj$Z2[, seq_len(ncol(obj$Z2) - (k - 1) * p) +
                                      (k - 1) * p, drop = F],
                             obj$Psi[, seq_len(ncol(obj$Z2) - (k - 1) * p) +
                                       (k - 1) * p, drop = F]
    )
  attributes(eps) <- attributes(res)
  ans <- rbind(init, rmfilter(eps, A, init))
  
  return(ans)
}

simulate.I2 <- function(obj, nsim, seed, res, init, ...) {
  stopifnot(!missing(res))
  
  p <- nrow(obj$alpha)
  r <- ncol(obj$alpha)
  k <- obj$lags
  
  A <- VAR(obj)
  
  if(missing(init)) init <- obj$X[seq_len(k), ]
  
  eps <- as.matrix(res) +
    tcrossprod(obj$Z2[, seq_len(ncol(obj$Z2) - p) + p, drop = F] %*%
                 obj$beta[seq_len(ncol(obj$Z2) - p) + p, , drop = F] +
                 obj$Z1[, seq_len(ncol(obj$Z1) - p) + p, drop = F] %*%
                 obj$delta[seq_len(ncol(obj$Z1) - p) + p, , drop = F],
               obj$alpha
    ) + obj$Z1[, seq_len(ncol(obj$Z1) - p) + p, drop = F] %*%
    tcrossprod(obj$tau[seq_len(ncol(obj$Z1) - p) + p,
                       , drop = F],
               obj$zeta
    )
  if (!is.null(ncol(obj$Z3)) && ncol(obj$Z3) > (k - 2) * p)
    eps <- eps + tcrossprod(obj$Z3[, seq_len(ncol(obj$Z3) - (k - 2) * p) +
                                     (k - 2) * p, drop = F],
                            obj$Psi[, seq_len(ncol(obj$Z3) - (k - 2) * p) +
                                      (k - 2) * p, drop = F]
    )
  attributes(eps) <- attributes(res)
  
  ans <- rbind(init, rmfilter(eps, A, init))
  
  return(ans)
}
