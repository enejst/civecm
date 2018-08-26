setrankI1 <- function(obj, r) {
  ## if (!inherits(obj, 'I1')) stop('Object must have class I1')
  ## if (!any(r %in% 0:ncol(obj$Z0))) stop('Illegal choice of rank')
  
  if (r > 0) {
    obj$beta <- obj$vectors[, seq_len(r), drop = FALSE]
    obj$alpha <- t(lm.fit(obj$R1 %*% obj$beta, obj$R0)$coef)
  }
  else {
    obj$alpha <- matrix(, ncol(obj$Z0), 0)
    obj$beta <- matrix(, ncol(obj$Z1), 0)
  }
  
  obj <- auxI1(obj)
  return(obj)
}

setrankI2 <- function(obj, r) {
  ## Z <- with(obj, merge(Z0, Z1, Z2, Z3))
  Time <- nrow(obj$Z0)
  
  if (!is.null(ncol(obj$Z3))) {
    R0 <- lm.fit(obj$Z3, obj$Z0)$residuals
    R1 <- lm.fit(obj$Z3, obj$Z1)$residuals
    R2 <- lm.fit(obj$Z3, obj$Z2)$residuals
  }
  else {
    R0 <- obj$Z0
    R1 <- obj$Z1
    R2 <- obj$Z2
  }
  
  ft <- unclass(obj)
  if (r[1] == 0) {
    if (r[2] < ncol(R0)) {
      tmpft <- rrr(R0, R1, , ncol(R0) - r[2])
      ft$tau <- tmpft$beta
      ft$rho <- matrix( , ncol(R0) - r[2], 0)
      ft$psi <- matrix( , ncol(R1), 0)
      ft$kappa <- lm.fit(R1 %*% ft$tau, R0)$coefficient
    } else {
      ft$tau <- matrix( , ncol(R1), 0)
      ft$rho <- matrix( , ncol(R0) - r[2], 0)
      ft$psi <- matrix( , ncol(R1), 0)
      ft$kappa <- matrix( , 0, ncol(R0))
    }
  } else {
    ## ft <- tauswitch(R0, R1, R2, r[1], r[2])
    p0 <- ncol(R0)
    p1 <- ncol(R1)
    s1 <- p0 - sum(r)
    ## initest <- rrr(R1[,1:p0], R2, NULL, r[1])
    initest <- rrr(diff(obj$X)[-seq_len(obj$lags),], obj$Z2, obj$Z1, r[1])
    ## initest <- rrr(diff(obj$X)[-seq_len(obj$lags),], lag(obj$X)[-seq(obj$lags),], , r[1])
    tAu <- initest$beta[, 1:r[1], drop = FALSE]
    if (s1 > 0) tAu	<- cbind(tAu, Null(tAu)[, 1:s1])
    tmpft <- .Fortran('tauswitch'
                      , as.double(R0)
                      , as.double(R1)
                      , as.double(R2)
                      , as.double(diag(p1))
                      , time = as.integer(Time)
                      , as.integer(p0)
                      , as.integer(p1)
                      , as.integer(r[1])
                      , as.integer(r[2])
                      , as.integer(p1)
                      , tau = tAu
                      , rho = matrix(0.0, r[1] + s1, r[1])
                      , psi = matrix(0.0, p1, r[1])
                      , kappa = matrix(0.0, r[1] + s1, p0 - r[1])
                      , as.integer(10000)
                      , 1e-9
                      , integer(1)
    )
    ft <- within(ft, {
      tau <- tmpft$tau
      rho <- tmpft$rho
      psi <- tmpft$psi
      kappa <- tmpft$kappa
      iter <- tmpft$time
    }
    )
    if (ft$iter == 0) stop('No convergence')
    if (tmpft[[17]] != 0) stop('Error in tauswitch')
  }
  
  ft$R0 <- R0
  ft$R1 <- R1
  ft$R2 <- R2
  ft <- auxI2(ft, r[1])
  
  ## ft$frequency <- obj$frequency
  return(ft)
}
