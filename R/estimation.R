rrr <- function(Z0, Z1, Z2, r) {
  Z0.loc <- as.matrix(Z0)
  Z1.loc <- as.matrix(Z1)
  if (missing(Z2))
    Z2.loc <- Z2 <- NULL
  else if (is.null(Z2))
    Z2.loc <- NULL
  else if (is.null(ncol(Z2)))
    Z2.loc <- NULL
  else
    Z2.loc <- as.matrix(Z2)
  
  Z <- cbind(Z0.loc, Z1.loc, Z2.loc)
  Time <- nrow(Z)
  M <- crossprod(Z) / Time
  if (is.null(Z2.loc)) {
    R0 <- Z0
    R1 <- Z1
  }
  else {
    R0 <- lm.fit(Z2, Z0)$residuals
    R1 <- lm.fit(Z2, Z1)$residuals
  }
  
  ## The following is to be implemented in Fortran
  if(!is.matrix(R0)) {R0 <- xts(as.matrix(coredata(R0),ncol=1),index(R0))}
  svdR0 <- svd(R0)
  svdR1 <- svd(R1)
  svd01 <- svd(crossprod(svdR0$u, svdR1$u))
  
  ans	        <- list()
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
  ans$R0      <- R0
  ans$R1      <- R1
  ans$beta    <- ans$vectors[, 1:r, drop = FALSE]
  ans$alpha   <- t(lm.fit(ans$R1%*%ans$beta, ans$R0)$coef)
  ans$Z0      <- Z0
  ans$Z1      <- Z1
  ans$Z2      <- Z2
  ans$M       <- M
  return(ans)
}


tauswitch <- function(R0, R1, R2, r, s2, Hmat = NULL) {
  R0 <- as.matrix(R0)
  R1 <- as.matrix(R1)
  R2 <- as.matrix(R2)
  
  Time <- nrow(R0)
  p <- ncol(R0)
  p1 <- ncol(R1)
  s1 <- p - r - s2
  
  if (!is.null(Hmat)) {
    if (nrow(Hmat) != p1 | ncol(Hmat) > p1) stop('Restriction matrix has wrong dimensions')
  }
  else Hmat <- diag(p1)
  
  initest <- rrr(R0, R1%*%Hmat, NULL, r + s1)
  tAu <- Hmat%*%initest$beta
  ## if(s1 > 0) tAu <- cbind(tAu, Null(tAu)[, 1:s1])
  
  ll <- -1e9
  while(1) {
    tmp.rrr	<- geigen(crossprod(lm.fit(cbind(R2 %*% tAu, R1), R0)$res) / Time,
                      crossprod(lm.fit(R1 %*% tAu, R0)$res) /Time)
    aLpha <- crossprod(lm.fit(R1 %*% tAu, R0)$res) %*% tmp.rrr$vectors[, (p - r + 1):p, drop = FALSE] / Time
    aLpha.ort <- tmp.rrr$vectors[, 0:(p - r), drop = FALSE]
    aLpha.bar <- t(solve(crossprod(aLpha), t(aLpha)))
    
    tmp1.r <- lm.fit(cbind(R0 %*% aLpha.ort, R2 %*% tAu, R1), R0 %*% aLpha.bar)
    rHo <- as.matrix(tmp1.r$coef)[(p - r + 1):(p + s1), , drop = FALSE]
    pSi <- as.matrix(tmp1.r$coef)[(p + s1 + 1):(p + p1 + s1), , drop = FALSE]
    ## NB! Division by Time us omited for simplicity
    OMega1 <- crossprod(tmp1.r$residuals)
    
    tmp2.r <- lm.fit(R1 %*% tAu, R0 %*% aLpha.ort)
    kAppa <- tmp2.r$coef
    OMega2 <- crossprod(tmp2.r$residuals)
    
    mA <- rHo %*% solve(OMega1, t(rHo))
    mB <- crossprod(lm.fit(cbind(R0 %*% aLpha.ort, R1), R2)$res%*%Hmat)
    mC <- kAppa  %*% solve(OMega2, t(kAppa))
    mD <- crossprod(R1%*%Hmat)
    mE <- (rHo %*% solve(OMega1,
                         crossprod(aLpha.bar,
                                   crossprod(lm.fit(cbind(R0 %*% aLpha.ort,
                                                          R1),
                                                    R0)$res,
                                             lm.fit(cbind(R0 %*% aLpha.ort,
                                                          R1), R2)$res)
                         )
    ) +
      kAppa %*% solve(OMega2,
                      crossprod(aLpha.ort,
                                crossprod(R0,
                                          R1)
                      )
      )
    ) %*% Hmat
    tAu0 <- tAu
    tAu <- Hmat %*% t(matrix(solve(mB %x% mA + mD %x% mC,
                                   as.numeric(mE)
    ),
    ncol = ncol(Hmat)
    )
    )
    
    ll0 <- ll
    ll <- -0.5 * Time * c(determinant(crossprod(lm.fit(cbind(R2 %*% tAu %*% rHo +
                                                               R1 %*% pSi,
                                                             R1%*%tAu),
                                                       R0)$res) / Time
    )$modulus
    )
    ## cat(ll,'\t')
    ## if(ll < ll0 - 0.0001) stop('Likelihood is decreasing')
    if (abs((ll - ll0) / ll0) < 1e-12) break
  }
  aLpha.ort.bar <- aLpha.ort%*%solve(crossprod(aLpha.ort))
  oMega <- t(as.matrix(tmp1.r)[0:(p - r), , drop = FALSE])
  return(list(alpha = aLpha, tau = tAu, rho = rHo, psi = pSi, kappa = kAppa))
}

I1gls <- function(obj, r) {
  ## Lütkepol's two step procedure for estimating the cointegrting relations
  cl <- match.call()
  tmp.fit <- lm.fit(obj$R1, obj$R0)
  Pi0	<- t(tmp.fit$coef)
  obj$Sigma <- crossprod(tmp.fit$residuals) / (nrow(obj$Z1) - ncol(obj$Z1))
  
  obj$alpha <- Pi0[, 1:r, drop = FALSE]
  
  bEta0 <- solve(crossprod(obj$R1[, -(1:r), drop = FALSE]),
                 crossprod(crossprod(obj$R0 -
                                       tcrossprod(obj$R1[, 1:r, drop = FALSE],
                                                  obj$alpha),
                                     obj$R1[, -(1:r), drop = FALSE]),
                           solve(obj$Sigma, obj$alpha))) %*%
    solve(crossprod(obj$alpha,
                    solve(obj$Sigma,
                          obj$alpha)
    )
    )
  
  obj$call <- cl
  obj$beta <- rbind(diag(1, r), bEta0)
  
  obj	<- auxI1(obj)
  class(obj) <- 'I1gls'
  
  return(obj)
}