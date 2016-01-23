
tracetest <- function(fit, type, a.p.value, bootstrap, reps, cluster) UseMethod('tracetest')

tracetest.I1 <- function(fit, type = 'I1', a.p.value = TRUE, bootstrap = FALSE, bootstrapType = "Resampling", reps = 499, cluster = NULL) {
  if (type == 'I1') {
    p0 <- ncol(fit$Z0)
    ## Johansen 1995, second edition
    ## MacKinnon et al 1998. JAE. UDEN TILLADELSE...
    ans <- data.frame(p.r = p0:1)
    ans$r <- 0:(p0 - 1)
    ans$eig.value <- fit$values[1:p0]
    ans$Trace <- -nrow(fit$Z0) * rev(cumsum(rev(log(1 - fit$values[1:p0]))))
    if (a.p.value) {
      if (!fit$call[["dettrend"]] %in% c('none', 'mean', 'drift')) stop('Deterministic specification not covered')
      ans	<- cbind(ans,
                   tracemckinnon(ans$p.r,
                                 'trace',
                                 fit$call[['dettrend']],
                                 c(0.05, 0.01)
                   ),
                   tracemckinnon(ans$p.r,
                                 'trace',
                                 fit$call[['dettrend']],
                                 tst = ans$Trace)
      )
    }
    if (bootstrap) {
      tmp0 <- numeric(nrow(ans))
      for (i in seq(nrow(ans))) {
        cat('Bootstrapping for rank', i - 1, 'of', nrow(ans) - 1, '\n')
        tmp1 <- unlist(boot(function(x) tracetest(x,
                                                  type = 'I1',
                                                  FALSE,
                                                  FALSE
        )$Trace[i],
        update(fit, i - 1),
        reps,
        cluster = cluster
        )
        )
        tmp0[i]	<- mean(tmp1 > ans$Trace[i])
      }
      ans	<- cbind(ans, p.value = tmp0)
    }
    ## Note: right now only stadard cases are covered with cirtical values and p values. Later bootstrap methods will be supplied.
    if (!is.null(fit$call[['exogenous']]))
      warning('Critical values are unreliable because of exogenous variables')
    ##    ans <- as.data.frame(ans)
    class(ans) <- 'tracetest.I1'
  }
  else if (type == 'I2') {
    p0 <- ncol(fit$Z0)
    if (p0 > 8)
      cat('Critical values only avalible for eigth strochastic trends\n')
    ans <- list()
    ans$values <- matrix(NA, p0, p0 + 1)
    ans$p.values <- ans$values
    HpOmega <- crossprod(update(fit, p0)$residuals) / nrow(fit$Z0)
    load(paste(system.file(package='civecm'), 'extdata/I2crit.Rdata', sep = '/'))
    for (i in seq_len(p0)) {
      for (j in i:(p0 + 1)) {
        tmpfit <- update(fit, c(i - 1, p0 - j + 1))
        H0Omega	<- crossprod(tmpfit$residuals) / nrow(fit$Z0)
        ans$values[i, j] <- -nrow(fit$Z0) * log(det(solve(H0Omega, HpOmega)))
        if (bootstrap) {
          cat('Bootstrapping rank', sprintf('(%d,%d)', i - 1, p0 - j + 1), '\n')
          tmp <- unlist(boot(function(x) -2 * (logLik(x) - logLik(update(x, p0))),
                             tmpfit,
                             reps,
                             cluster = cluster))
          ans$p.values[i, j] <- mean(ans$values[i, j] < tmp)
        } else if (p0 < 9) {
          tmp <- c(I2crit[I2crit$s1 == (j - i) & I2crit$s2 == (p0 - j + 1),9:10])
          ans$p.values[i, j] <- pgamma(ans$values[i, j],
                                       shape = tmp$mean^2 / tmp$var,
                                       scale = tmp$var / tmp$mean,
                                       lower.tail = F
          )
        }
      }
    }
    class(ans) <- 'tracetest.I2'
  }
  else stop('Type must be either I1 or I2')
  return(ans)
}

tracetest2.I1 <- function(fit, type = 'I1', a.p.value = TRUE, boottype = "Resampling", reps = 499, cluster = NULL) {
  if (type == 'I1') {
    p0 <- ncol(fit$Z0)
    ## Johansen 1995, second edition
    ## MacKinnon et al 1998. JAE. UDEN TILLADELSE...
    ans <- data.frame(p.r = p0:1)
    ans$r <- 0:(p0 - 1)
    ans$eig.value <- fit$values[1:p0]
    ans$Trace <- -nrow(fit$Z0) * rev(cumsum(rev(log(1 - fit$values[1:p0]))))
    if (a.p.value) {
      if (!fit$call[["dettrend"]] %in% c('none', 'mean', 'drift')) stop('Deterministic specification not covered')
      ans	<- cbind(ans,
                   tracemckinnon(ans$p.r,
                                 'trace',
                                 fit$call[['dettrend']],
                                 c(0.05, 0.01)
                   ),
                   tracemckinnon(ans$p.r,
                                 'trace',
                                 fit$call[['dettrend']],
                                 tst = ans$Trace)
      )
    }
    if (bootstrap) {
      tmp0 <- numeric(nrow(ans))
      for (i in seq(nrow(ans))) {
        cat('Bootstrapping for rank', i - 1, 'of', nrow(ans) - 1, '\n')
        tmp1 <- unlist(boot(function(x) tracetest(x,
                                                  type = 'I1',
                                                  FALSE,
                                                  FALSE
        )$Trace[i],
        update(fit, i - 1),
        reps,
        cluster = cluster
        )
        )
        tmp0[i]	<- mean(tmp1 > ans$Trace[i])
      }
      ans	<- cbind(ans, p.value = tmp0)
    }
    ## Note: right now only stadard cases are covered with cirtical values and p values. Later bootstrap methods will be supplied.
    if (!is.null(fit$call[['exogenous']]))
      warning('Critical values are unreliable because of exogenous variables')
    ##    ans <- as.data.frame(ans)
    class(ans) <- 'tracetest.I1'
  }
  else if (type == 'I2') {
    p0 <- ncol(fit$Z0)
    if (p0 > 8)
      cat('Critical values only avalible for eigth strochastic trends\n')
    ans <- list()
    ans$values <- matrix(NA, p0, p0 + 1)
    ans$p.values <- ans$values
    HpOmega <- crossprod(update(fit, p0)$residuals) / nrow(fit$Z0)
    load(paste(system.file(package='civecm'), 'extdata/I2crit.Rdata', sep = '/'))
    for (i in seq_len(p0)) {
      for (j in i:(p0 + 1)) {
        tmpfit <- update(fit, c(i - 1, p0 - j + 1))
        H0Omega	<- crossprod(tmpfit$residuals) / nrow(fit$Z0)
        ans$values[i, j] <- -nrow(fit$Z0) * log(det(solve(H0Omega, HpOmega)))
        if (bootstrap) {
          cat('Bootstrapping rank', sprintf('(%d,%d)', i - 1, p0 - j + 1), '\n')
          tmp <- unlist(boot(function(x) -2 * (logLik(x) - logLik(update(x, p0))),
                             tmpfit,
                             reps,
                             cluster = cluster))
          ans$p.values[i, j] <- mean(ans$values[i, j] < tmp)
        } else if (p0 < 9) {
          tmp <- c(I2crit[I2crit$s1 == (j - i) & I2crit$s2 == (p0 - j + 1),9:10])
          ans$p.values[i, j] <- pgamma(ans$values[i, j],
                                       shape = tmp$mean^2 / tmp$var,
                                       scale = tmp$var / tmp$mean,
                                       lower.tail = F
          )
        }
      }
    }
    class(ans) <- 'tracetest.I2'
  }
  else stop('Type must be either I1 or I2')
  return(ans)
}

tracetest.I2 <- function(fit, type = 'I2', a.p.value = TRUE, bootstrap = FALSE, reps = 100, cluster = NULL) NextMethod(update(fit, ncol(fit$X)))

tracemckinnon <- function(npr, type = c('trace', 'lambda'), dettrend = c('none', 'mean', 'drift'), tlev, tst) {
  ## Function is just an interface to McKinnons Fortran function. Remember to ask
  ## for accept of use.
  
  npr <- if (is.vector(npr) & all(npr %in% 0:12)) as.integer(npr)
  else stop('npr must be a vector of integer values between 0 and 12')
  
  itt <- if (type == 'lambda') as.integer(1)
  else if (type == 'trace') as.integer(2)
  else stop('Test type must be either trace or lambda')
  
  itv <- if (dettrend == 'none') as.integer(0)
  else if (dettrend == 'mean') as.integer(3)
  else if (dettrend == 'drift') as.integer(4)
  else stop('Deterministic specification must be either none, mean or trend')
  
  if (!missing(tlev)) {
    arg <- if (all(0 < tlev) & all(tlev < 0.5)) as.numeric(tlev)
    else stop('Testlevel must be between 0.0001 and 0.5')
    nc <- as.integer(1)
  }
  else if (!missing(tst)) {
    arg <- if (all(0 < tst) & length(tst) == length(npr)) as.numeric(tst)
    else stop('All teststatistic must be positive and in same number as restrictions')
    nc <- as.integer(2)
  }
  val <- numeric(1)
  rpath <- system.file(package='civecm')
  
  if (nc == 1) {
    ans <- matrix(NA, length(npr), length(arg))
    colnames(ans) <- paste('Cr', substr(arg, 2, 5), sep = '')
    for (i in seq_along(npr)) {
      for (j in seq_along(arg)) {
        ans[i, j] <- .Fortran('johval',
                              npr[i],
                              itt,
                              itv,
                              nc,
                              isave = integer(1),
                              arg[j],
                              val = val,
                              rpath = rpath,
                              as.integer(nchar(rpath))
        )$val
      }
    }
  }
  else if (nc == 2) {
    ans <- matrix(NA, length(npr), 1)
    colnames(ans) <- 'a.p.value'
    for (i in seq_along(npr)) {
      ans[i, 1] <- .Fortran('johval',
                            npr[i],
                            itt,
                            itv,
                            nc,
                            isave = integer(1),
                            arg[i],
                            val = val,
                            rpath = rpath,
                            as.integer(nchar(rpath))
      )$val
    }
  }
  return(ans)
}
