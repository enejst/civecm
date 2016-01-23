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

testAutocor <- function(obj, portmanteau = TRUE, lagrange = c(1, 2)) {
  sm <- summary(obj)
  if (!is.null(obj$H.alpha) || !is.null(obj$Ha))
    warning('Portmanteau test is invalid when alpha is restricted.')
  
  Time <- nrow(obj$Z0)
  test <- character(0)
  
  ## Portmanteau test, as in Lütkepohl (2007). Originally from Hisking (1980)
  if (!portmanteau) {
    lb.value <- numeric(0)
    lb.df    <- numeric(0)
  }
  else {
    test     <- c(test, 'Portmanteau')
    lb.value <- 0
    h        <- floor(Time / 4)
    for (i in seq_len(h)) {
      OMega.h  <- crossprod(residuals(obj)[(1 + i):Time, ],
                            resid(obj)[1:(Time - i), ]) / Time
      lb.value <- lb.value + sum(diag(solve(sm$Omega, t(OMega.h)) %*%
                                        solve(sm$Omega, OMega.h))) / (Time - i)
    }
    
    lb.value <- lb.value * Time^2
    lb.df <- sm$p['p0']^2 * h -
      sm$p['p0']^2 * (sm$lags - 1) -
      sm$p['p0'] * ncol(obj$beta)
  }
  
  ## LM test, as in Lütkepohl (2007)
  lm.value <- numeric(0)
  lm.df <- numeric(0)
  
  if (any(lagrange > 0)) {
    lagrange <- sort(lagrange)
    test <- c(test, sprintf('LM(%d)', lagrange))
    V <- with(obj, merge(xts(Z1 %*% beta, index(Z1)), Z2))
    
    for (i in seq_len(max(lagrange))) {
      V <- merge(V, lag(sm$residuals, i))
      if (i %in% lagrange) {
        V[is.na(V)] <- 0
        tmp.res <- sm$residuals - xts(V%*%solve(crossprod(V),
                                                crossprod(V, sm$residuals)
        ), index(V))
        OMega.e	<- crossprod(tmp.res) / Time
        lm.value <- c(lm.value, Time * (sm$p['p0'] - sum(diag(solve(sm$Omega,
                                                                    OMega.e
        )
        )
        )
        )
        )
        lm.df <- c(lm.df, i * sm$p['p0']^2)
      }
    }
  }
  
  df <- c(lb.df, lm.df)
  value <- c(lb.value, lm.value)
  
  ans <- data.frame(Type = test,
                    Distr = 'ChiSq',
                    Df = df,
                    Value = value,
                    p.value = 1 - pchisq(value, df)
  )
  return(ans)
}

testExclude <- function(obj) {
  if (!inherits(obj, 'I1')) stop('Object must have class I1')
  
  r.max <- ncol(obj$X) - 1
  values <- p.vals <- matrix(NA, r.max, nrow(obj$beta))
  tmp <- matrix('', 2 * r.max, nrow(obj$beta) + 2)
  
  for (r in seq_len(r.max)) {
    fb <- update(obj, r)
    for (i in seq_len(nrow(obj$beta))) {
      H <- diag(nrow(obj$beta))[, -i]
      tmp <- anova.I1(fb, restrictBeta(fb, H), df = r)
      values[r, i] <- tmp[1, 1]
      p.vals[r, i] <- tmp[1, 3]
    }
  }
  colnames(values) <- colnames(obj$Z1)
  out <- list(df = 1:r.max, value = values, p.value = p.vals)
  class(out) <- 'autotest.I1'
  
  return(out)
}

testLagLength <- function(fit, lag.max = 8) {
  
  if (class(fit) != 'I1') stop('Input must have class I1')
  
  ##	m <- match(c('data', 'lags'), names(fit$call))
  Time <- nrow(fit$X) - lag.max
  
  loglik <- Regr <- rep(NA, lag.max)
  infocr <- matrix(NA, lag.max, 2)
  
  for (i in lag.max:1) {
    fit$call[['data']] <- if (i == lag.max) fit$X else fit$X[-seq_len(lag.max - i),]
    fit$call[['lags']] <- i
    tmp <- eval(fit$call, envir = attr(fit$call, 'environment'))
    loglik[i] <- logLik(tmp)
    Regr[i] <- ncol(tmp$Z1) + ifelse(is.null(tmp$Z2), 0, ncol(tmp$Z2))
    infocr[i,] <- infocrit(tmp)[1:2]
  }
  
  z <- list()
  
  z$model.sum	<- data.frame(model = paste('VAR(', lag.max:1, ')', sep = ''),
                            k	= lag.max:1,
                            t	= Time,
                            Regr = rev(Regr),
                            loglik = formatC(rev(loglik), digits = 3, format = 'f'),
                            SC = formatC(rev(infocr[,1]), digits = 3, format = 'f'),
                            HQ = formatC(rev(infocr[,2]), digits = 3, format = 'f')
  )
  
  lagtest <- data.frame()
  for (i in lag.max:1) {
    for (j in (lag.max - 1):1) {
      if (i > j) {
        test <- sprintf('VAR(%i) << VAR(%i)', j, i)
        typetmp	<- ncol(tmp$X) * (Regr[i] - Regr[j])
        type <- sprintf('ChiSqr(%i)', typetmp)
        value <- 2 * (loglik[i] - loglik[j])
        pvalue <- 1 - pchisq(value, typetmp)
        
        z$lagtest <- rbind(z$lagtest,
                           data.frame(test,
                                      type,
                                      value = formatC(value, digits = 3, format = 'f'),
                                      p.value = formatC(pvalue, digits = 3, format = 'f')
                           )
        )
      }
    }
  }
  
  ##    return(lapply(z, format, digits = 6, scientific = F))
  return(z)
}

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

testResiduals <- function(obj) {
  sm <- summary(obj)
  ans <- list()
  
  ans$cor <- sm$Omega / sqrt(diag(sm$Omega) %o% diag(sm$Omega))
  ans$se <- sqrt(diag(sm$Omega))
  
  ans$info <- infocrit(obj)
  ans$autocor <- testAutocor(obj)
  ans$normal <- testNormal(residuals(obj))
  ans$arch <- testARCH(obj)
  
  return(ans)
}

testStationary <- function(obj, incl.trend = TRUE, incl.exo = TRUE, incl.sd = TRUE) {
  if(!inherits(obj, 'I1')) stop('Object must have class I1')
  
  H.index  <- ncol(obj$X)
  det.text <- NULL
  p1       <- ncol(obj$Z1)
  
  if ('exogenous' %in% names(obj$call) && incl.exo) {
    H.index	 <- c(H.index, 1:ncol(as.matrix(obj$exogenous)) +
                    H.index[length(H.index)])
    det.text <- c(det.text, 'Exogenous variables included in test\n')
    p1       <- p1 - ncol(as.matrix(obj$exogenous))
  }
  
  if ('shift.dummy' %in% names(obj$call) && incl.sd) {
    H.index	 <- c(H.index, 1:ncol(as.matrix(obj$shift.dummy)) +
                    H.index[length(H.index)])
    det.text <- c(det.text, 'Shift dummies included in test\n')
    p1       <- p1 - ncol(as.matrix(obj$shift.dummy))
  }
  
  if (obj$call$dettrend != 'none' && incl.trend) {
    H.index	 <- c(H.index, 1 + H.index[length(H.index)])
    det.text <- c(det.text, 'Trend included in test\n')
    p1       <- p1 - 1
  }
  
  r.max <- ncol(obj$X) - 1
  
  values <- p.vals <- matrix(NA, r.max, ncol(obj$X))
  tmp	<- matrix('', 2 * r.max, ncol(obj$X) + 2)
  
  for (r in seq_len(r.max)) {
    fb <- eval(obj$call, envir = attr(obj$call, 'environment'))
    fb <- update(fb, r)
    
    for (i in seq_len(ncol(obj$X)))	{
      H.index[1]   <- i
      H            <- lapply(1:r, function(x) diag(nrow(obj$beta))[,-i])
      H[[1]]       <- diag(nrow(obj$beta))[, H.index, drop = F]
      tmp          <- anova.I1(fb, restrictBeta(fb, H), df = p1 - r)
      values[r, i] <- tmp[1, 1]
      p.vals[r, i] <- tmp[1, 3]
    }
  }
  
  colnames(values) <- colnames(obj$X)
  out	             <- list(df = ncol(obj$Z1) - 1:r.max, value = values, p.value = p.vals)
  class(out)       <- 'autotest.I1'
  return(out)
}

testUnit <- function(obj) {
  if(!inherits(obj, 'I1')) stop('Object must have I1')
  
  r.max <- ncol(obj$X) - 1
  
  values <- p.vals <- matrix(NA, r.max, ncol(obj$X))
  tmp <- matrix('', 2 * r.max, nrow(obj$alpha))
  
  for (r in seq_len(r.max)) {
    fb <- update(obj, r)
    for (i in seq_len(nrow(obj$alpha))) {
      H <- diag(nrow(obj$alpha))[, i, drop = FALSE]
      tmp <- anova.I1(fb,
                      restrictAlpha(fb, H, ifelse(r != 1, 'known', 'restr')),
                      df = ncol(obj$X) - r
      )
      values[r, i] <- tmp[1, 1]
      p.vals[r, i] <- tmp[1, 3]
    }
  }
  
  colnames(values) <- colnames(obj$X)
  out <- list(df = ncol(obj$X) - 1:r.max, value = values, p.value = p.vals)
  class(out) <- 'autotest.I1'
  return(out)
}

testWeakexo <- function(obj) {
  if(!inherits(obj, 'I1')) stop('Object must have class I1')
  
  r.max <- ncol(obj$X) - 1
  
  values <- p.vals <- matrix(NA, r.max, ncol(obj$X))
  tmp	<- matrix('', 2 * r.max, nrow(obj$alpha))
  
  for (r in seq_len(r.max)) {
    fb <- update(obj, r)
    for (i in seq_len(nrow(obj$alpha))) {
      H <- as.matrix(diag(nrow(obj$alpha))[, -i])
      tmp <- anova.I1(fb, restrictAlpha(fb, H), df = r * (ncol(obj$X) - ncol(H)))
      values[r, i] <- tmp[1, 1]
      p.vals[r, i] <- tmp[1, 3]
    }
  }
  
  colnames(values) <- colnames(obj$X)
  out <- list(df = 1:r.max, value = values, p.value = p.vals)
  class(out) <- 'autotest.I1'
  return(out)
}