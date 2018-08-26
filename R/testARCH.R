testARCH <- function(obj, order) {
  if (!inherits(obj, 'I1')) stop('Object must have class I1')
  
  i.k <- obj$call[['lags']]
  i.p <- ncol(obj$Z0)
  i.T <- nrow(obj$Z0)
  if (missing(order)) i.q <- i.k
  else i.q <- order
  
  m.res <- na.omit(lag(xts(t(apply(resid(obj), 1, tcrossprod)), index(obj$Z0)), seq(0,i.q)))
  m.Sigma <- crossprod(lm(m.res[,1:i.p^2]~m.res[,-(1:i.p^2)])$residuals) / i.T
  m.0Sigma <- crossprod(scale(m.res[,1:i.p^2], scale = F)) / i.T
  f.R2 <- 1 - 2 / i.p / (i.p + 1) * sum(diag(ginv(m.0Sigma) %*% m.Sigma))
  
  f.xts <- 0.5 * i.T * i.p * (i.p + 1) * f.R2
  d.xts <- data.frame(Distr = 'ChiSq',
                      Df = 0.25 * i.q * i.p^2 * (i.p + 1)^2,
                      Value = f.xts,
                      p.value = 1 - pchisq(f.xts, 0.25 * i.q * i.p^2 * (i.p + 1)^2)
  )
  
  f.uts <- rep(NA, i.p)
  for (i in seq_len(i.p)) {
    tmp <- na.omit(lag(obj$residuals[, i], 0:i.k))^2
    f.uts[i] <- summary(lm(tmp[, 1]~tmp[, -1]))$r.squared * (i.T - i.k)
  }
  
  d.uts <- data.frame(Variable = colnames(obj$residuals),
                      Distr = 'ChiSq',
                      Df = i.k,
                      Value = f.uts,
                      p.value = 1 - pchisq(f.uts, i.k)
  )
  
  return(list(Univariate.test = d.uts, Multivariate.test = d.xts))
}
