rrr <- function(Z0, Z1, Z2, r) {
  Z0.loc <- as.matrix(Z0)
  Z1.loc <- as.matrix(Z1)
  if (missing(Z2))
    Z2.loc <- Z2 <- NULL
  else if (is.null(Z2))
    Z2.loc <- NULL
  else
    Z2.loc <- as.matrix(Z2)
  Z <- cbind(Z0.loc, Z1.loc, Z2.loc)
  
  Time <- nrow(Z)
  p0 <- ncol(Z0.loc)
  p1 <- ncol(Z1.loc)
  #	p2		<- ncol(Z2.loc)
  
  M <- crossprod(Z) / Time
  
  if (is.null(Z2.loc)) {
    R0 <- Z0
    R1 <- Z1
  }
  else {
    R0 <- lm.fit(Z2, Z0)$residuals
    R1 <- lm.fit(Z2, Z1)$residuals
  }
  ##	ans <- .Fortran('rrr', as.double(R0), as.double(R1), as.integer(Time),
  ## as.integer(p0), as.integer(p1), values = double(p1),
  ## vectors = matrix(0.0, p1, p1))
  
  ## The following is to be implemented in Fortran
  svdR0 <- svd(R0)
  svdR1 <- svd(R1)
  svd01 <- svd(crossprod(svdR0$u, svdR1$u))
  ans	<- list()
  ans$values	<- svd01$d^2
  ans$vectors	<- svd01$v * sqrt(Time)
  for (i in 1:length(svdR1$d))
  {
    if (svdR1$d[i] > 2e-16)
      ans$vectors[i,] <- ans$vectors[i,]/svdR1$d[i]
    else
      ans$vectors[i,] <- 0
  }
  ans$vectors <- svdR1$v%*%ans$vectors
  ans$R0 <- R0
  ans$R1 <- R1
  ans$S11 <- crossprod(R1) / Time
  ans$beta <- ans$vectors[, 1:r, drop = FALSE]
  ans$alpha <- t(lm.fit(ans$R1%*%ans$beta, ans$R0)$coef)
  ans$Z0 <- Z0
  ans$Z1 <- Z1
  ans$Z2 <- Z2
  ans$M <- M
  return(ans)
}
