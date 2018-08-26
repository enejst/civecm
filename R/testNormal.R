
testNormal <- function(obj) {
  if (inherits(obj, 'I1')) res <- resid(obj)
  else res <- obj
  ## Doornik and Hansen 2008 Oxford Bulletin of Economics and Statistics
  
  i.T <- nrow(res)
  m.Omega <- crossprod(res) / i.T
  m.C <- m.Omega / sqrt(diag(m.Omega) %o% diag(m.Omega))
  l.eig.C <- eigen(m.C)
  m.H <- l.eig.C$vectors
  m.Lamb <- diag(l.eig.C$values)
  m.u <- scale(res, scale = FALSE) %*% diag(diag(m.Omega)^(-0.5)) %*% m.H %*% tcrossprod(diag(diag(m.Lamb)^(-0.5)), m.H)
  
  f.b1.sqr <- colMeans(m.u^3)
  f.b2 <- colMeans(m.u^4)
  
  f.delta <- (i.T - 3) * (i.T + 1) * (i.T^2 + 15 * i.T - 4) * 6
  f.a	<- (i.T - 2) * (i.T + 5) * (i.T + 7) * (i.T^2 + 27 * i.T - 70) / f.delta
  f.c	<- (i.T - 7) * (i.T + 5) * (i.T + 7) * (i.T^2 + 2 * i.T - 5) / f.delta
  f.k	<- (i.T + 5) * (i.T + 7) * (i.T^3 + 37 * i.T^2 + 11 * i.T - 313) * 0.5 / f.delta
  f.alpha <- f.a + f.b1.sqr^2 * f.c
  f.chi <- (f.b2 - 1 - f.b1.sqr^2) * 2 * f.k
  f.z2 <- ((0.5 * f.chi / f.alpha)^(1/3) - 1 + 1 / 9 / f.alpha) * 3 * sqrt(f.alpha)
  
  f.beta <- 3 * (i.T^2 + 27 * i.T - 70) * (i.T + 1) * (i.T + 3) / (i.T - 2) / (i.T + 5) / (i.T + 7) / (i.T + 9)
  f.omega2 <- -1 + sqrt(2 * (f.beta - 1))
  f.delta <- log(sqrt(f.omega2))^(-0.5)
  f.y <- f.b1.sqr * sqrt((f.omega2 - 1) * (i.T + 1) * (i.T + 3) / 12 / (i.T - 2))
  f.z1 <- f.delta * log(f.y + sqrt(f.y^2 + 1))
  
  f.unorm.test <- f.z1^2 + f.z2^2
  d.unorm.test <- data.frame(Variable = colnames(res),
                             Distr = 'ChiSq',
                             Df = 2,
                             Value = f.unorm.test,
                             p.value = 1 - pchisq(f.unorm.test, 2)
  )
  f.mnorm.test <- sum(f.unorm.test)
  d.mnorm.test <- data.frame(Distr = 'ChiSq',
                             Df = 2*ncol(m.u),
                             Value = f.mnorm.test,
                             p.value = 1 - pchisq(f.mnorm.test, 2 * ncol(m.u))
  )
  
  f.skew <- skewness(res)
  f.kurt <- kurtosis(res)
  
  d.uni.stat <- data.frame(Mean = colMeans(res),
                           Std.dev = sqrt(diag(m.Omega)),
                           Skewness = f.skew,
                           Kurtosis = f.kurt
  )
  
  ans	<- list(Univariate.test = d.unorm.test,
              Multivariate.test = d.mnorm.test,
              Summary.stat = d.uni.stat
  )
  
  return(ans)
}
