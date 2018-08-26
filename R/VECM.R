VECM <- function(data, lags, dettrend, rank, season, dummy, exogenous, method = 'rrr') {
  cl <- match.call()
  mf <- match.call()
  
  m <- match(c('data', 'dummy', 'exogenous'), names(mf), 0L)
  if (any(missing(data), missing(lags), missing(dettrend))) stop('The arguments data, lags and dettrend are required')
  ## if (!inherits(data, 'xts')) stop('Data must be of class xts')
  if (any(is.na(data))) stop('Data contains NAs')
  if (!any(lags %in% seq(12))) stop ('The model must have lags between 1 and 12')
  ## if (!missing(exogenous)) {
  ##     if (!inherits(exogenous, 'xts')) stop('Exogenous variables must be of class
  ## 	 xts')
  ##     if (sum(is.na(exogenous))) stop('Exogenous variables contains NAs')
  ## }
  ## if (!missing(dummy)) {
  ##     if (!inherits(dummy, 'xts')) stop('Dummy variables must be of class xts')
  ##     if (any(is.na(dummy))) stop('Dummy variables contains NAs')
  ## }
  ## if (!dettrend %in% c('none', 'mean', 'drift')) stop('The determistic specification must be either none, mean or drift')
  
  mf[[1L]] <- as.name("CIModelFrame")
  mf <- eval(mf, parent.frame())
  
  if (missing(rank)) {
    ft <- rrr(mf$Z0, mf$Z1, mf$Z2, ncol(mf$Z0))
    ft$lags <- lags
    ft <- auxI1(ft)
  }
  else if (length(rank) == 1) {
    ft <- rrr(mf$Z0, mf$Z1, mf$Z2, ncol(mf$Z0))
    ft$lags <- lags
    ft <- setrankI1(ft, rank)
  }
  else if (length(rank) == 2) {
    ## if (!(cl[['dettrend']] %in% c('none', 'drift'))) stop('Deterministic specification must me either none og drift')
    if (lags < 2) stop('There must be at least tow lags in the model')
    mf$lags <- lags
    ft <- setrankI2(mf, rank)
  }
  else stop('Illegal choice of rank')
  
  ft$call <- cl
  attr(ft$call, 'environment') <- parent.frame()
  ft$X <- mf$X
  
  return(ft)
}

