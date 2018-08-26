
plot.eigenCompanion <- function(x, y, ...) {
  stopifnot (inherits(x, 'eigenCompanion'))
  
  plot(x$values, asp = 1, xlim = c(-1, 1), ylim = c(-1, 1), xlab = 'Real', ylab = 'Imaginary')
  symbols(0, 0, circles = 1, inches = FALSE, add = TRUE)
  abline(h = 0, v = 0)
}

plot.I1	<- function(x, y, type = 'CI.rel', ...) {
  if (missing(y)) type = 'CI.rel'
  else if (missing(type)) type = y
  if (type == 'CI.rel') {
    Z1beta <- xts(x$Z1 %*% x$beta, index(x$Z1))
    R1beta <- xts(x$R1 %*% x$beta, index(x$R1))
    op <- par(ask = TRUE, mfrow = c(2, 1))
    for (i in seq_len(ncol(x$beta))) {
      plot(Z1beta[,i],
           ylab = paste("beta", i, "'Z1", sep = ''),
           main = paste('Cointegration relation', i))
      plot(R1beta[,i],
           ylab = paste("beta", i, "'R1", sep = ''),
           main = '')
    }
    par(op)
  }
  else if (type == 'residuals') {
    op <- par(ask = TRUE, mfrow = c(2, 2))
    for (i in seq_len(ncol(x$X))) {
      plot(x$Z0[, i],
           col = 2,
           ylab = '',
           main = paste('Actual and fitted of', colnames(x$Z0)[i]))
      lines(x$fitted.values[, i],
            col = 4)
      acf(ts(x$residuals[, i]),
          main = paste('Residual autocorrelogram of', colnames(x$Z0)[i]))
      plot(rstandard(x)[, i],
           ylab = '',
           main = paste('Standardized residuals of', colnames(x$Z0)[i]))
      plot(density(rstandard(x)[, i]),
           ylab = '',
           col = 2,
           main = paste('Estimated and standard normal\n density of',
                        colnames(x$Z0)[i],
                        'residuals'))
      curve(dnorm,
            add = TRUE,
            col = 4)
    }
    par(op)
  }
  else stop("Type must be either CI.rel or residuals")
}

plot.I2 <- function(x, y, type, ...) {
  if (missing(y)) type = 'CI.rel'
  else if (missing(type)) type = y
  if (type == 'CI.rel') {
    ZCIm <- with(x, xts(Z2 %*% beta + Z1 %*% delta, index(Z2)))
    RCIm <- with(x, xts(R2 %*% beta + R1 %*% delta, index(R2)))
    op <- par(ask = TRUE, mfrow = c(2,1))
    for (i in seq_len(ncol(ZCIm))) {
      plot(ZCIm[,i],
           ylab = paste("beta", i, "'Z2 + delta", i, "'Z1", sep = ''),
           main = paste('Multi-cointegration relation', i))
      plot(RCIm[,i],
           ylab = paste("beta", i, "'R2 + delta", i, "'R1", sep = ''),
           main = '')
    }
    par(op)
    NextMethod()
  }
  else NextMethod()
}
