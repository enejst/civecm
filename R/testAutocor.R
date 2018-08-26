testAutocor <- function(obj, portmanteau = TRUE, lagrange = c(1, 2)) {
  sm <- summary(obj)
  if (!is.null(obj$H.alpha))
    warning('Portmanteau test is invalid then alpha is restricted.')
  Time <- nrow(obj$Z0)
  
  test <- character(0)
  
  ## Portmanteau test, as in Lütkepohl (2007). Originally from Hisking (1980)
  if (!portmanteau) {
    lb.value <- numeric(0)
    lb.df <- numeric(0)
  }
  else {
    test <- c(test, 'Portmanteau')
    lb.value <- 0
    h <- floor(Time / 4)
    for (i in seq_len(h)) {
      OMega.h <- crossprod(residuals(obj)[(1 + i):Time, ],
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
