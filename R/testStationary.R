testStationary <- function(obj, incl.trend = TRUE, incl.exo = TRUE, incl.sd = TRUE) {
  if(!inherits(obj, 'I1')) stop('Object must have class I1')
  
  H.index <- ncol(obj$X)
  det.text <- NULL
  p1 <- ncol(obj$Z1)
  
  if ('exogenous' %in% names(obj$call) && incl.exo) {
    H.index	<- c(H.index, 1:ncol(as.matrix(obj$exogenous)) +
                   H.index[length(H.index)])
    det.text <- c(det.text, 'Exogenous variables included in test\n')
    p1 <- p1 - ncol(as.matrix(obj$exogenous))
  }
  
  if ('shift.dummy' %in% names(obj$call) && incl.sd) {
    H.index	<- c(H.index, 1:ncol(as.matrix(obj$shift.dummy)) +
                   H.index[length(H.index)])
    det.text <- c(det.text, 'Shift dummies included in test\n')
    p1 <- p1 - ncol(as.matrix(obj$shift.dummy))
  }
  
  if (obj$call$dettrend != 'none' && incl.trend) {
    H.index	<- c(H.index, 1 + H.index[length(H.index)])
    det.text <- c(det.text, 'Trend included in test\n')
    p1 <- p1 - 1
  }
  
  r.max <- ncol(obj$X) - 1
  
  values <- p.vals <- matrix(NA, r.max, ncol(obj$X))
  tmp	<- matrix('', 2 * r.max, ncol(obj$X) + 2)
  
  for (r in seq_len(r.max)) {
    fb <- eval(obj$call, envir = attr(obj$call, 'environment'))
    fb <- update(fb, r)
    
    for (i in seq_len(ncol(obj$X)))	{
      H.index[1] <- i
      H <- lapply(1:r, function(x) diag(nrow(obj$beta))[,-i])
      H[[1]] <- diag(nrow(obj$beta))[, H.index, drop = F]
      tmp <- anova.I1(fb, restrictBeta(fb, H), df = p1 - r)
      values[r, i] <- tmp[1, 1]
      p.vals[r, i] <- tmp[1, 3]
    }
  }
  
  colnames(values) <- colnames(obj$X)
  out	<- list(df = ncol(obj$Z1) - 1:r.max, value = values, p.value = p.vals)
  class(out) <- 'autotest.I1'
  return(out)
}
