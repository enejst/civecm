CIModelFrame <- function(data, lags, dettrend, rank, season = FALSE, exogenous, dummy, ...) {
  X <- try.xts(data)
  if (is.null(colnames(X))) colnames(X) <- paste('X', seq(ncol(X)), sep='.')
  k <- lags
  fullSampleVec <- index(X)
  fullSample <- paste(format(as.POSIXct(first(fullSampleVec)), '%Y-%m-%d'),
                      format(as.POSIXct(last(fullSampleVec)), '%Y-%m-%d'),
                      sep = '/')
  effectiveSample <- paste(format(as.POSIXct(first(fullSampleVec[-seq_len(k)])), '%Y-%m-%d'),
                           format(as.POSIXct(last(fullSampleVec)), '%Y-%m-%d'),
                           sep = '/')
  Time <- nrow(X) - k
  if (dettrend == 0) {
    xtsTime <- xts(c(rep(0, k), rep(1, Time)), fullSampleVec)
    colnames(xtsTime) <- 'constant'
  }
  else if(dettrend > 0) {
    xtsTime <- xts(cbind(1, poly(seq_len(Time + k), dettrend, raw = T)), fullSampleVec)
    colnames(xtsTime) <- c('constant', 'linear', 'quadratic', 'cubic', ifelse(dettrend > 3, NULL, paste('poly', seq(4, dettrend), sep = '.')))[seq_len(dettrend + 1)]
  }
  
  if (missing(exogenous)) W <- l.d.W <- l.d2.W <- NULL
  else {
    W <- try.xts(exogenous)
    if (is.null(colnames(W))) colnames(W) <- paste('exo', seq_len(ncol(W)), sep = '.')
  }
  
  
  if (!season) Season <- NULL
  else {
    periodX <- periodicity(X)$scale
    if (periodX == 'monthly') freq <- 12
    else if (periodX == 'quarterly') freq <- 4
    else stop('Seasonal factors are not supplied for other periodicities than month and quarterly')
    Season <- xts(matrix(rep(c(1, rep(0, freq - 2), -1), Time + k - 1)[seq_len((freq - 1) * (Time + k))], Time + k, byrow = T), fullSampleVec)
    colnames(Season) <- paste('Season', seq_len(ncol(Season)), sep = '.')
  }
  
  if (missing(dummy)) Dummy <- NULL
  else Dummy <- try.xts(dummy)
  
  d.X	<- diff(X)
  colnames(d.X) <- paste('d', colnames(X), sep = '.')
  l.X	<- lag(X)
  colnames(l.X) <- colnames(X)
  
  if (missing(rank) || length(rank) == 1) {
    if(k == 1) l.d.X <- xts(, fullSampleVec)
    else {
      l.d.X <- lag(d.X, seq(k - 1))
      colnames(l.d.X) <- paste('l',
                               rep(seq(k - 1),
                                   each = ncol(X)),
                               '.',
                               colnames(d.X), sep = '')
    }
    
    if (is.null(ncol(W))) l.d.W <- NULL
    else {
      l.d.W <- lag(diff(W), seq_len(k) - 1)
      colnames(l.d.W) <- c(paste('d', colnames(W), sep = '.'),
                           if (k == 1) NULL
                           else paste('l', rep(seq(k - 1), each = ncol(W)), '.d.', colnames(W), sep = ''))
    }
    
    if (dettrend == -1)
      Det1 <- Det2 <- NULL
    else if (dettrend == 0) {
      Det1 <- xtsTime[, 1]
      Det2 <- NULL
    }
    else {
      Det1 <- xtsTime[, dettrend + 1]
      Det2 <- xtsTime[, seq_len(dettrend)]
    }
    
    tmp <- list(X = X)
    tmp$Z0 <- d.X[effectiveSample]
    tmp$Z1 <- merge(l.X, W, Det1)[effectiveSample]
    Z2 <- merge(l.d.X, l.d.W, Season, Dummy, Det2)
    if (is.null(ncol(Z2))) tmp$Z2 <- l.d.X
    else tmp$Z2 <- Z2[effectiveSample]
    stopifnot(all(sapply(tmp, function(x) all(!is.na(x)))))
    class(tmp) <- 'model.frame.I1'
  }
  else if (length(rank) == 2) {
    l.d.X <- lag(diff(X))
    colnames(l.d.X) <- paste('l', '.d.', colnames(X), sep = '')
    d2.X <- diff(d.X)
    colnames(d2.X) <- paste('d2', colnames(X), sep = '.')
    if (k == 1) stop('In the I(2) model there must be at least two lags')
    else if (k == 2) l.d2.X <- xts(, fullSampleVec)
    else {
      l.d2.X <- lag(d2.X, seq_len(k - 2))
      colnames(l.d2.X) <- paste('l', rep(1:(k - 2), each = ncol(X)), '.d2.', colnames(X), sep = '')
    }
    
    if (is.null(ncol(W))) d.W <- l.d2.W <- NULL
    else {
      d.W <- diff(W)
      l.d2.W <- lag(diff(d.W), seq_len(k - 1) - 1)
      colnames(d.W) <- paste('d', colnames(W), sep = '.')
      colnames(l.d2.W) <- c(paste('d2', colnames(W), sep = '.'),
                            if (k == 2) NULL
                            else paste('l', rep(seq_len(k - 2), each = ncol(W)), '.d2.', colnames(W), sep =''))
    }
    
    if (dettrend == -1)
      Det1 <- Det2 <- Det3 <- NULL
    else if (dettrend == 0) {
      Det2 <- xtsTime[, 1]
      Det1 <- Det3 <- NULL
    }
    else if (dettrend == 1) {
      Det2 <- xtsTime[, 2]
      Det1 <- xtsTime[, 1]
      Det3 <- NULL
    }
    else {
      Det2 <- xtsTime[, dettrend + 1]
      Det1 <- xtsTime[, dettrend]
      Det3 <- xtsTime[, seq_len(dettrend - 1)]
    }
    
    tmp <- list(X = X)
    tmp$Z0 <- d2.X[effectiveSample]
    tmp$Z1 <- merge(l.d.X, d.W, Det1)[effectiveSample]
    tmp$Z2 <- merge(l.X, W, Det2)[effectiveSample]
    Z3 <- merge(l.d2.X, l.d2.W, Season, Dummy, Det3)
    if (is.null(ncol(Z3))) tmp$Z3 <- l.d2.X
    else tmp$Z3 <- Z3[effectiveSample]
    stopifnot(all(sapply(tmp, function(x) all(!is.na(x)))))
    class(tmp) <- 'model.frame.I2'
  }
  
  return(tmp)
}
