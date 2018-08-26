print.autotest.I1 <- function(x, digits = 3, ...) {
  out	<- coefmatrix(x$value, x$p.value, digits = digits)
  lout <- matrix(, 0, 3)
  colnames(lout) <- c('r', 'df', '5pct.c.v.')
  for (i in seq(length(x$df))) {
    lout <- rbind(lout, c(i, x$df[i], formatC(qchisq(0.95, x$df[i]),
                                              digits = digits,
                                              format = 'f'
    )
    ), ""
    )
  }
  out	<- cbind(lout, out)
  print(out, right = TRUE, quote = FALSE)
}

print.eigenCompanion <- function(x, ...) {
  stopifnot (inherits(x, 'eigenCompanion'))
  eigFrame <- data.frame(Real = Re(x$values),
                         Imaginary = Im(x$values),
                         Modulus = Mod(x$values),
                         Argument = Arg(x$values)
  )
  print(eigFrame, digits = 3)
}

print.I1 <- function(x, ...) {
  cat('\nCall:\n')
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  
  if (ncol(x$beta) == 0)
    cat('\nModel has rank zero. No long run parameters to print\n')
  else {
    cat('\nbeta (transposed)\n')
    print(t(x$beta))
  }
}

print.I1gls <- function(x, ...) {
  cat('\nCall:\n')
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  
  if (is.null(x$beta))
    cat('\nModel has rank zero. No long run parameters to print\n')
  else {
    cat('\nBeta (transposed)\n')
    print(t(x$beta))
  }
}

print.I2 <- function(x, ...) {
  cat('\nCall:\n')
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  
  ## if(ncol(x$beta) == 0) cat('\nModel has rank zero. No long run parameters to print\n')
  ## else
  ## {
  cat('\nbeta (transposed)\n')
  print(t(x$beta))
  cat('\nbeta1 (transposed)\n')
  print(t(x$beta1))
  ## }
}

print.summary.I1 <- function(x, ...) {
  cat('\nCall:\n')
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  
  periodX <- periodicity(x$X)
  periodZ <- periodicity(x$Z0)
  cat('\nEffective Sample: ',
      paste(periodZ$start, collapse = ':'),
      ' to ',
      paste(periodZ$end, collapse = ':'),
      ' (',
      nrow(x$residuals),
      ' observations)',
      sep = ''
  )
  cat('\nFrequency:       ', periodZ$scale, '\n')
  
  if (ncol(x$beta) == 0)
    cat('\nModel has rank zero. No long run parameters to print\n')
  else {
    cat('\nbeta (transposed)\n')
    print(t(x$beta))
    
    cat('\nalpha:\n')
    print(coefmatrix(x$alpha, x$t.alpha), right = TRUE, quote = FALSE)
    
    cat('\nPi:\n')
    print(coefmatrix(x$Pi, x$t.Pi), right = TRUE, quote = FALSE)
  }
  cat('\nLog-Likelihood: ', x$logLik, '\n')
}

print.summary.I2 <- function(x, ...) {
  cat('Nothing yet, so here is the normal output\n')
  print.I2(x)
}

print.tracetest.I1 <- function(x, ...) {
  cat('\nThe I(1) Tracetest Statistics:\n')
  ans <- as.data.frame(unclass(x))
  ans[, -seq_len(2)] <- lapply(ans[, -seq_len(2)], function(y) formatC(y, digits = 3, format = 'f'))
  print(ans, row.names = FALSE)
}

print.tracetest.I2 <- function(x, ...) {
  ans	<- coefmatrix(x$values, x$p.values, sig.dig = F)
  ans <- sub('\\(NA\\)', '    ', ans)
  ans	<- sub('NA', '  ', ans)
  ans	<- cbind(c(rbind(seq(nrow(x$values)) - 1, ' ')), ans)
  ans	<- rbind(c(' ', 's2', rep('', ncol(x$values) - 1)),
               c('r', rev(seq(ncol(x$values)) - 1)), ans)
  cat('\nThe I(2) Tracetest Statistics:')
  prmatrix(ans, right = TRUE, quote = FALSE, collab = rep('', ncol(ans)),
           rowlab = rep('', nrow(ans))
  )
}
