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
