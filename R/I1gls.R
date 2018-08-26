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
