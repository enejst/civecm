#
# An R package for cointegration
# Author: Andreas Noack Jensen
# Date: April 2010
#
# Notes:
# The tracetest p-values has been stolen from McKinnon and I need to ask him
# in case I need to redistribute this code.
#  - We will simulate our own!
#
# ---------------------------------------------
# Date 09/12/2015
# Development taken over by Emil Nejstgaard
#

anova.I1 <- function(objList, ..., df) {
  # Construct a Likelihood Ratio test between two models. The test
  # assumes a chi-square limiting distribution and is hence not 
  # adapted to nonstandard cases.
  #
  # Args
  #   objList: a list of objects of type I1
  #   df: the number of degrees of freedom the test has
  # 
  # Return:
  #   data.frame containing the test value, the coresponding p-value from 
  #   a chi-sqaure distribution with df degrees of freedom
  #
  
  if(length(objects <- list(objList, ...)) == 1) stop('Nothing to test against')
  test <- 2 * (logLik(objects[[1]]) - logLik(objects[[2]]))
  return(data.frame(Value = test,
                    Test = paste('ChiSq(', df, ')', sep = ''),
                    p.value = 1 - pchisq(test, df)
  )
  )
}

CIModelFrame <- function(data, lags, dettrend, rank, season = FALSE, exogenous, ...) {
  # Function that creates a model.frame object. Either I1 or I2. 
  # 
  # Args:
  #   data: an xts object with data for analysis
  #   lags: number of lags in the model
  #   dettrend: specified the deterministic trend structure
  #       - none
  #       - mean
  #       - drift
  #       - qtrend
  #   rank: cointegration rank
  #   season: Boolean for seasonals, month and quarter provided
  #   exogenous: an xts object with exogenous variables (Optional)
  #       (Note - exogenous varibles are treatet as standard variables in
  #       with regard to lagging and differencing. Only they are not
  #       included on the left hand side in the VAR equation.)
  #       Dummy variables should be included through the exogenous variables.
  #     
  # Return:
  #   model.frame Object
  #     - model.frame.I1 (central object for I1 models)
  #     - model.frame.I2 (central object for I2 models)
  #
  
  k               <- lags
  X               <- try.xts(data)
  if (is.null(colnames(X))) colnames(X) <- paste('X', seq(ncol(X)), sep='.')
  
  fullSampleVec   <- index(X)
  effectiveSample <- fullSampleVec[-seq_len(k)]
  Time <- nrow(X) - k
  
  # Check if exogenous variables have been provided
  if (missing(exogenous) || is.null(exogenous)) W <- l.W <- l.d.W <- l.d2.W <- NULL
  else {
    W <- try.xts(exogenous)
    if (is.null(colnames(W))) colnames(W) <- paste('exo', seq_len(ncol(W)), sep = '.')
    
    if(k == 1) {
      l.W           <- lag(W)
      colnames(l.W) <- colnames(W)
      l.d.W         <- NULL
    }
    else {
      l.W             <- lag(W)
      colnames(l.W)   <- colnames(W)
      d.W             <- diff(W,1)
      l.d.W           <- lag(d.W, seq(k - 1))
      colnames(l.d.W) <- paste('l',
                               rep(seq(k - 1), each = ncol(W)),
                               '.',
                               colnames(d.W),
                               sep = '')
    }
  }
  
  # Construct deterministic trends
  trendList <- dettrendbuild(fullSampleVec,lags,dettrend)
  l.detTrend <- trendList[["l.detTrend"]]
  l.d.detTrend <- trendList[["l.d.detTrend"]]
  l.d2.detTrend <- trendList[["l.d2.detTrend"]]
  
  # Check for inclusion of seasonals
  if (!season) Season <- NULL
  else {
    periodX <- periodicity(X)$scale
    if (periodX == 'monthly') freq <- 12
    else if (periodX == 'quarterly') freq <- 4
    else stop('Seasonal factors are not supplied for other periodicities than month and quarterly')
    Season <- xts(matrix(rep(c(1, 
                               rep(0, freq - 2), 
                               -1), 
                             Time + k - 1)[seq_len((freq - 1) * (Time + k))], 
                         Time + k, byrow = T), 
                  fullSampleVec)
    
    colnames(Season) <- paste('Season', seq_len(ncol(Season)), sep = '.')
  }
  
  # Constructing lags and first differences
  d.X	          <- diff(X)
  colnames(d.X) <- paste('d', colnames(X), sep = '.')
  l.X	          <- lag(X)
  colnames(l.X) <- paste('l', colnames(l.X), sep='.')
  
  if (missing(rank) || length(rank) == 1) {
    # Include lagged differences, NA if number of lags does not exceed k
    if(k == 1) 
      l.d.X <- NULL
    else {
      l.d.X           <- lag(d.X, seq(k - 1))
      colnames(l.d.X) <- paste('l',
                               rep(seq(k - 1), each = ncol(X)),
                               '.',
                               colnames(d.X),
                               sep = '')
    }
    
    
    tmp    <- list(X = X)
    tmp$Z0 <- d.X[effectiveSample]
    tmp$Z1 <- merge(l.X, l.W, l.detTrend)[effectiveSample]
    Z2 <- merge(l.d.X, l.d.W, l.d.detTrend, Season)
    if (is.null(ncol(Z2)) || ncol(Z2)==0) tmp$Z2 <- l.d.X
    else tmp$Z2 <- Z2[effectiveSample]
    stopifnot(all(sapply(tmp, function(x) all(!is.na(x)))))
    class(tmp) <- 'model.frame.I1'
  }
  else if (length(rank) == 2) {
    l.d.X           <- lag(diff(X))
    colnames(l.d.X) <- paste('l', '.d.', colnames(X), sep = '')
    d2.X            <- diff(d.X)
    colnames(d2.X)  <- paste('d2', colnames(X), sep = '.')
    if (k == 1) stop('In the I(2) model there must be at least two lags')
    else if (k == 2) l.d2.X <- xts(, fullSampleVec)
    else {
      l.d2.X <- lag(d2.X, seq_len(k - 2))
      colnames(l.d2.X) <- paste('l', rep(1:(k - 2), each = ncol(X)), '.d2.', colnames(X), sep = '')
    }
    
    if (is.null(ncol(W))) l.W <- l.d.W <- l.d2.W <- NULL
    else {
      l.W   <- lag(W)
      colnames(l.W) <- colnames(W)
      l.d.W <- diff(l.W)
      colnames(l.d.W) <- paste('d', colnames(W), sep = '.')
      l.d2.W <- lag(diff(d.W), seq_len(k - 1) - 1)
      colnames(d.W) <- paste('d', colnames(W), sep = '.')
      colnames(l.d2.W) <- c(paste('d2', colnames(W), sep = '.'),
                            if (k == 2) NULL
                            else paste('l', rep(seq_len(k - 2), each = ncol(W)), '.d2.', colnames(W), sep =''))
    }
    
    tmp    <- list(X = X)
    tmp$Z0 <- d2.X[effectiveSample]
    tmp$Z1 <- merge(l.d.X, l.d.W, l.d.detTrend)[effectiveSample]
    tmp$Z2 <- merge(l.X, l.W, l.detTrend)[effectiveSample]
    Z3     <- merge(l.d2.X, l.d2.W, l.d2.detTrend, Season)
    if (is.null(ncol(Z3))) tmp$Z3 <- l.d2.X
    else tmp$Z3 <- Z3[effectiveSample]
    stopifnot(all(sapply(tmp, function(x) all(!is.na(x)))))
    class(tmp) <- 'model.frame.I2'
  }
  
  return(tmp)
}

coefmatrix <- function(coef, stat1, digits = 3, sig.dig = c(TRUE, FALSE)) {
  # Print the coefficient matrices along with standar errors and t
  # statistics
  
  m.coef  <- as.matrix(coef)
  m.stat1 <- as.matrix(stat1)
  mchar   <- max(nchar(cbind(m.coef, m.stat1)))
  if (any(dim(m.coef)!=dim(m.stat1))) stop('Coef and stat matrices must have same dimension')
  if (length(sig.dig) != 2)
    if (length(sig.dig) == 1) sig.dig <- rep(sig.dig, 2)
  else stop('Argument sig.dig must have length two')
  m.out <- matrix(, 0, ncol(coef))
  for (i in seq_len(nrow(coef))) {
    m.out <- rbind(m.out,
                   formatC(m.coef[i, ],
                           width = -1,
                           digits = digits,
                           format = ifelse(sig.dig[1], 'fg', 'f')
                   ),
                   paste('(',
                         formatC(m.stat1[i, ],
                                 width = -1,
                                 digits = digits,
                                 format = ifelse(sig.dig[2], 'fg', 'f')),
                         ')',
                         sep = '')
    )
  }
  maxi    <- apply(m.out, 2, nchar)
  li      <- apply(m.out, 2, function(x) regexpr('.', x, fixed = TRUE))
  ri      <- maxi - li
  maxcolr <- apply(ri, 2, max)
  maxcoll <- apply(li, 2, max)
  for (i in seq_len(nrow(m.out))) {
    for (j in seq_len(ncol(m.out))) {
      m.out[i, j]	<- paste(paste(rep(' ', maxcoll[j] - li[i, j]),
                                 collapse = ''
      ),
      m.out[i, j],
      paste(rep(' ', maxcolr[j] - ri[i, j]),
            collapse = ''
      ),
      sep = ''
      )
    }
  }
  if (!is.null(rownames(coef))) rownames(m.out) <- c(rbind(rownames(coef), ""))
  colnames(m.out) <- colnames(coef)
  return(m.out)
}

companion <- function(object) {
  stopifnot(inherits(object, 'VAR'))
  
  if (dim(object)[3] > 1)
    ans <- rbind(matrix(c(object), dim(object)[1], dim(object)[1] * dim(object)[3]),
                 cbind(diag((dim(object)[3] - 1) * dim(object)[1]),
                       matrix(0, (dim(object)[3] - 1) * dim(object)[1],
                              ifelse(dim(object)[3] > 1, dim(object)[1], 0))))
  else ans <- object[,,1]
  return(ans)
}

datestring <- function(dates, freq) {
  ans <- if(freq %in% c(4, 12)) sapply(dates, function(x) paste(floor(x), ifelse(freq == 4, ' Q', ' M'), round((x - floor(x))*freq + 1), sep = ''))
  else dates
  return(ans)
}

dettrendbuild <- function(dateIndex,lags,dettrend) {
  # Construct the deterministic trend to include in the system
  # 
  # Args:
  #   dateIndex: an index of the dates matrich the dates in the
  #               effective sample
  #   
  #   lags: number of lags in the model
  #
  #   dettrend: specification of the trend 
  #       - none 
  #       - mean
  #       - drift
  #       - qtrend
  #
  # Returns
  #   trendList: list with the deterministic trends to included
  #              inside and outside the cointegration relations
  # 
  
  # Construct deterministic trends
  if (dettrend == "mean") {
    detTrend             <- xts(rep(1, length(dateIndex)), dateIndex)
    l.detTrend           <- lag(detTrend)[-seq_len(lags)]
    colnames(l.detTrend) <- 'ci-constant'
    l.d.detTrend         <- NULL
    l.d2.detTrend        <- NULL
  }
  else if(dettrend == "urmean") {
    detTrend             <- xts(rep(1, length(dateIndex)), dateIndex)
    l.detTrend           <- NULL
    l.d.detTrend         <- lag(detTrend)[-seq_len(lags)]
    colnames(l.detTrend) <- 'constant'
    l.d2.detTrend        <- NULL
  }
  else if(dettrend == "drift") {
    # Construct xts variable with a constant and polynomial trends
    detTrend               <- xts(cbind(poly(seq_len(length(dateIndex)), 1, raw = T)), dateIndex)
    l.detTrend             <- lag(detTrend)
    colnames(l.detTrend)   <- c('ci-trend')
    if(lags>1) {
      l.d.detTrend <- diff(lag(detTrend))
      colnames(l.d.detTrend) <- "constant"
    }
    else {
      l.d.detTrend <- NULL
    }
    l.d2.detTrend          <- NULL
  } 
  else if(dettrend == "qtrend" && lags > 1) {
    detTrend                <- xts(cbind(poly(seq_len(length(dateIndex)), 2, raw = T)), dateIndex)
    l.detTrend              <- lag(detTrend)
    colnames(l.detTrend)    <- c('ci-qtrend')
    l.d.detTrend            <- diff(detTrend)
    colnames(l.d.detTrend)  <- c('trend')
    l.d2.detTrend           <- diff(diff(detTrend))
    colnames(l.d2.detTrend) <- c('constant')
  }
  else {
    l.detTrend <- NULL
    l.d.detTrend <- NULL
    l.d2.detTrend <- NULL
  }
  
  trendList <- list('l.detTrend'=l.detTrend,
                    'l.d.detTrend'=l.d.detTrend,
                    'l.d2.detTrend'=l.d2.detTrend)
}

hessian <- function(obj) {
  # Calculates the covariance matrix for the parameter
  # estimates given restrictions imposed through the function
  # restrictLongRun
  #
  # Args:
  #   obj: an object of type I1
  #
  # Return
  #   covarList: a list with covariance matrices and 
  #
  
  p0     <- ncol(obj$Z0)
  p1     <- ncol(obj$Z1)
  p2     <- ifelse(is.null(obj$Z2), NULL, ncol(obj$Z2))
  G      <- if(is.null(obj$Ha)) { diag(p0*r) } else {obj$Ha}    
  H      <- if(is.null(obj$Hb)) { diag(p1*r) } else {obj$Hb}
  alpha  <- obj$alpha
  beta   <- obj$beta
  Omega  <- obj$Omega
  iOmega <- solve(Omega)
  r      <- ncol(alpha)
  T      <- nrow(obj$Z1)
  S11    <- crossprod(obj$R1,obj$R1) / T
  
  hessian_11 <- t(G) %*% kronecker(iOmega, crossprod(beta,S11) %*% beta) %*% G
  hessian_12 <- t(G) %*% kronecker(iOmega %*% alpha, crossprod(beta,S11)) %*% H 
  hessian_22 <- t(H) %*% kronecker(crossprod(alpha,iOmega) %*% alpha, S11) %*% H
  
  hessian <- T * rbind(cbind(hessian_11,hessian_12),
                       cbind(t(hessian_12),hessian_22))
  iHessian <- solve(hessian)
  
  sdPsi   <- if(ncol(G)==1) {
    sqrt(iHessian[1:ncol(G),1:ncol(G)])
  }
  else {
    sqrt(diag(iHessian[1:ncol(G),1:ncol(G)]))}
  sdPhi   <- if(ncol(H)==1) {
    sqrt(iHessian[(ncol(G)+1):(ncol(G)+ncol(H)),(ncol(G)+1):(ncol(G)+ncol(H))])} 
  else {
    sqrt(diag(iHessian[(ncol(G)+1):(ncol(G)+ncol(H)),(ncol(G)+1):(ncol(G)+ncol(H))]))
  }
  
  sdAlpha <- matrix(t(G %*% sdPsi),nrow=p0,ncol=r)
  sdBeta  <- matrix(H %*% sdPhi,nrow=p1,ncol=r)
  vSdPi   <- kronecker(beta,diag(p0)) %*% vec(sdAlpha)
  sdPi    <- abs(matrix(vSdPi,nrow=p0,ncol=p1))
  
  covarList           <- list(hessian=hessian)
  covarList$iHessian  <- solve(iHessian)
  covarList$sePi      <- sdPi / sqrt(T)
  covarList$sdAlpha   <- sdAlpha
  covarList$sdBeta    <- sdBeta
  covarList$seAlpha   <- sdAlpha / sqrt(T)
  covarList$seBeta    <- sdBeta / sqrt(T)
  
  return(covarList)
}

I1gls <- function(obj, r) {
  ## Lütkepol's two step procedure for estimating the cointegrting relations
  cl <- match.call()
  tmp.fit <- lm.fit(obj$R1, obj$R0)
  Pi0	<- t(tmp.fit$coef)
  obj$Sigma <- crossprod(tmp.fit$residuals) / (nrow(obj$Z1) - ncol(obj$Z1))
  
  obj$alpha <- Pi0[, 1:r, drop = FALSE]
  
  bEta0 <- solve(crossprod(obj$R1[, -(1:r), drop = FALSE]),
                 crossprod(crossprod(obj$R0 -
                                       tcrossprod(obj$R1[, 1:r, drop = FALSE],
                                                  obj$alpha),
                                     obj$R1[, -(1:r), drop = FALSE]),
                           solve(obj$Sigma, obj$alpha))) %*%
    solve(crossprod(obj$alpha,
                    solve(obj$Sigma,
                          obj$alpha)
    )
    )
  
  obj$call <- cl
  obj$beta <- rbind(diag(1, r), bEta0)
  
  obj	<- auxI1(obj)
  class(obj) <- 'I1gls'
  
  return(obj)
}

infocrit <- function(obj) {
  if (!inherits(obj, 'I1')) stop('Object must be of typee either I1 or summary.I1')
  
  sm    <- summary(obj)
  Time  <- nrow(sm$residuals)
  logOm <- as.double(determinant(sm$Omega)$modulus)
  sc    <- as.double(logOm + (sm$p['p2'] * sm$p['p0'] +
                                sm$rank * (sm$p['p0'] + sm$p['p1'] - sm$rank)
  ) * log(Time) / Time)
  hq <- as.double(logOm + (sm$p['p2'] * sm$p['p0'] +
                             sm$rank * (sm$p['p0'] + sm$p['p1'] - sm$rank)
  ) * 2 * log(log(Time)) / Time)
  tc <- as.double(1 - sum(diag(solve(var(obj$Z0), sm$Omega))) / sm$p['p0'])
  
  return(c(SC = sc, HQ = hq, TC = tc))
}

largeres <- function(fit, abs.value) {
  res <- rstandard(fit)
  
  if (missing(abs.value)) idx <- cbind(apply(abs(res), 2, which.max), seq(ncol(res)))
  else idx <- which(abs(res) > abs.value, T)
  dates <- index(res)[idx[,1]]
  variable <- apply(idx, 1, function(x) colnames(res)[x[2]])
  values <- as.matrix(res)[idx]
  
  return(data.frame(dates = dates, var = variable, values = values))
}

logLik.I1 <- function(object, const = TRUE, ...) {
  OMega <- with(object, crossprod(residuals) / nrow(Z0))
  ll <- -nrow(object$Z0) * 0.5 * c(determinant(OMega)$modulus)
  if (const) ll <- ll - prod(dim(object$Z0)) * 0.5 * (log(2 * pi) + 1)
  names(ll) <- 'logLik'
  return(ll)
}

Null <- function(M) {
  # Function borrowed from MASS but with here is the convention that the Null of a zero column matrix
  # is the identity
  tmp <- qr(M)
  set <- if (tmp$rank == 0) seq_len(ncol(M))
  else -seq_len(tmp$rank)
  ans <- if (ncol(M) == 0) diag(nrow(M))
  else qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
  return(ans)
}

residuals.I1	<- function(object, ...) {
  tmp	<- object$residuals
  return(tmp)
}

rmfilter <- function(x, coef, init) {
  if (any(is.na(x))) stop('Series includes NAs')
  if (any(is.na(coef))) stop('Coefficients includes NAs')
  if (is.null(dim(x))) stop('Input series has no dimensions')
  nl <- dim(x)
  if (!is.array(coef)) stop('coef must be an array')
  cdim <- dim(coef)
  if (nl[2] != cdim[2] | cdim[1] != cdim[2]) stop('coef has wrong dimensions')
  if (missing(init)) init <- matrix(0, cdim[3], nl[2])
  else if (any(is.na(init))) stop('Initial values includes NAs')
  else if (any(dim(init) != c(cdim[3], nl[2]))) stop('init matrix has wrong dimensions')
  storage.mode(x) <- 'double'
  ans	<- .Fortran('rmfilter',
                  x,
                  as.integer(nl),
                  as.double(coef),
                  as.integer(cdim),
                  as.double(init)
  )[[1]]
  return(ans)
}

rrr <- function(Z0, Z1, Z2, r) {
  Z0.loc <- as.matrix(Z0)
  Z1.loc <- as.matrix(Z1)
  if (missing(Z2))
    Z2.loc <- Z2 <- NULL
  else if (is.null(Z2))
    Z2.loc <- NULL
  else if (is.null(ncol(Z2)))
    Z2.loc <- NULL
  else
    Z2.loc <- as.matrix(Z2)
  
  Z <- cbind(Z0.loc, Z1.loc, Z2.loc)
  Time <- nrow(Z)
  M <- crossprod(Z) / Time
  if (is.null(Z2.loc)) {
    R0 <- Z0
    R1 <- Z1
  }
  else {
    R0 <- lm.fit(Z2, Z0)$residuals
    R1 <- lm.fit(Z2, Z1)$residuals
  }
  
  ## The following is to be implemented in Fortran
  if(!is.matrix(R0)) {R0 <- xts(as.matrix(coredata(R0),ncol=1),index(R0))}
  svdR0 <- svd(R0)
  svdR1 <- svd(R1)
  svd01 <- svd(crossprod(svdR0$u, svdR1$u))
  
  ans	        <- list()
  ans$values	<- svd01$d^2
  ans$vectors	<- svd01$v * sqrt(Time)
  for (i in 1:length(svdR1$d))
  {
    if (svdR1$d[i] > 2e-16)
      ans$vectors[i,] <- ans$vectors[i,]/svdR1$d[i]
    else
      ans$vectors[i,] <- 0
  }
  ans$vectors <- svdR1$v%*%ans$vectors
  ans$R0      <- R0
  ans$R1      <- R1
  ans$beta    <- ans$vectors[, 1:r, drop = FALSE]
  ans$alpha   <- t(lm.fit(ans$R1%*%ans$beta, ans$R0)$coef)
  ans$Z0      <- Z0
  ans$Z1      <- Z1
  ans$Z2      <- Z2
  ans$M       <- M
  return(ans)
}

rstandard.I1 <- function(model, ...) {
  return(model$residuals / matrix(sqrt(diag(summary(model)$Omega)),
                                  nrow(model$residuals),
                                  ncol(model$residuals),
                                  byrow = T
  )
  )
}

setrankI1 <- function(obj, r) {
  ## if (!inherits(obj, 'I1')) stop('Object must have class I1')
  ## if (!any(r %in% 0:ncol(obj$Z0))) stop('Illegal choice of rank')
  
  if (r > 0) {
    obj$beta <- obj$vectors[, seq_len(r), drop = FALSE]
    obj$alpha <- t(lm.fit(obj$R1 %*% obj$beta, obj$R0)$coef)
  }
  else {
    obj$alpha <- matrix(, ncol(obj$Z0), 0)
    obj$beta <- matrix(, ncol(obj$Z1), 0)
  }
  
  obj <- auxI1(obj)
  return(obj)
}

setrankI2 <- function(obj, r) {
  ## Z <- with(obj, merge(Z0, Z1, Z2, Z3))
  Time <- nrow(obj$Z0)
  
  if (!is.null(ncol(obj$Z3))) {
    R0 <- lm.fit(obj$Z3, obj$Z0)$residuals
    R1 <- lm.fit(obj$Z3, obj$Z1)$residuals
    R2 <- lm.fit(obj$Z3, obj$Z2)$residuals
  }
  else {
    R0 <- obj$Z0
    R1 <- obj$Z1
    R2 <- obj$Z2
  }
  
  ft <- unclass(obj)
  if (r[1] == 0) {
    if (r[2] < ncol(R0)) {
      tmpft <- rrr(R0, R1, , ncol(R0) - r[2])
      ft$tau <- tmpft$beta
      ft$rho <- matrix( , ncol(R0) - r[2], 0)
      ft$psi <- matrix( , ncol(R1), 0)
      ft$kappa <- lm.fit(R1 %*% ft$tau, R0)$coefficient
    } else {
      ft$tau <- matrix( , ncol(R1), 0)
      ft$rho <- matrix( , ncol(R0) - r[2], 0)
      ft$psi <- matrix( , ncol(R1), 0)
      ft$kappa <- matrix( , 0, ncol(R0))
    }
  } else {
    ## ft <- tauswitch(R0, R1, R2, r[1], r[2])
    p0 <- ncol(R0)
    p1 <- ncol(R1)
    s1 <- p0 - sum(r)
    ## initest <- rrr(R1[,1:p0], R2, NULL, r[1])
    initest <- rrr(diff(obj$X)[-seq_len(obj$lags),], obj$Z2, obj$Z1, r[1])
    ## initest <- rrr(diff(obj$X)[-seq_len(obj$lags),], lag(obj$X)[-seq(obj$lags),], , r[1])
    tAu <- initest$beta[, 1:r[1], drop = FALSE]
    if (s1 > 0) tAu	<- cbind(tAu, Null(tAu)[, 1:s1])
    tmpft <- .Fortran('tauswitch'
                      , as.double(R0)
                      , as.double(R1)
                      , as.double(R2)
                      , as.double(diag(p1))
                      , time = as.integer(Time)
                      , as.integer(p0)
                      , as.integer(p1)
                      , as.integer(r[1])
                      , as.integer(r[2])
                      , as.integer(p1)
                      , tau = tAu
                      , rho = matrix(0.0, r[1] + s1, r[1])
                      , psi = matrix(0.0, p1, r[1])
                      , kappa = matrix(0.0, r[1] + s1, p0 - r[1])
                      , as.integer(10000)
                      , 1e-9
                      , integer(1)
    )
    ft <- within(ft, {
      tau <- tmpft$tau
      rho <- tmpft$rho
      psi <- tmpft$psi
      kappa <- tmpft$kappa
      iter <- tmpft$time
    }
    )
    if (ft$iter == 0) stop('No convergence')
    if (tmpft[[17]] != 0) stop('Error in tauswitch')
  }
  
  ft$R0 <- R0
  ft$R1 <- R1
  ft$R2 <- R2
  ft <- aux.I2(ft, r[1])
  
  ## ft$frequency <- obj$frequency
  return(ft)
}

tauswitch <- function(R0, R1, R2, r, s2, Hmat = NULL) {
  R0 <- as.matrix(R0)
  R1 <- as.matrix(R1)
  R2 <- as.matrix(R2)
  
  Time <- nrow(R0)
  p <- ncol(R0)
  p1 <- ncol(R1)
  s1 <- p - r - s2
  
  if (!is.null(Hmat)) {
    if (nrow(Hmat) != p1 | ncol(Hmat) > p1) stop('Restriction matrix has wrong dimensions')
  }
  else Hmat <- diag(p1)
  
  initest <- rrr(R0, R1%*%Hmat, NULL, r + s1)
  tAu <- Hmat%*%initest$beta
  ## if(s1 > 0) tAu <- cbind(tAu, Null(tAu)[, 1:s1])
  
  ll <- -1e9
  while(1) {
    tmp.rrr	<- geigen(crossprod(lm.fit(cbind(R2 %*% tAu, R1), R0)$res) / Time,
                      crossprod(lm.fit(R1 %*% tAu, R0)$res) /Time)
    aLpha <- crossprod(lm.fit(R1 %*% tAu, R0)$res) %*% tmp.rrr$vectors[, (p - r + 1):p, drop = FALSE] / Time
    aLpha.ort <- tmp.rrr$vectors[, 0:(p - r), drop = FALSE]
    aLpha.bar <- t(solve(crossprod(aLpha), t(aLpha)))
    
    tmp1.r <- lm.fit(cbind(R0 %*% aLpha.ort, R2 %*% tAu, R1), R0 %*% aLpha.bar)
    rHo <- as.matrix(tmp1.r$coef)[(p - r + 1):(p + s1), , drop = FALSE]
    pSi <- as.matrix(tmp1.r$coef)[(p + s1 + 1):(p + p1 + s1), , drop = FALSE]
    ## NB! Division by Time us omited for simplicity
    OMega1 <- crossprod(tmp1.r$residuals)
    
    tmp2.r <- lm.fit(R1 %*% tAu, R0 %*% aLpha.ort)
    kAppa <- tmp2.r$coef
    OMega2 <- crossprod(tmp2.r$residuals)
    
    mA <- rHo %*% solve(OMega1, t(rHo))
    mB <- crossprod(lm.fit(cbind(R0 %*% aLpha.ort, R1), R2)$res%*%Hmat)
    mC <- kAppa  %*% solve(OMega2, t(kAppa))
    mD <- crossprod(R1%*%Hmat)
    mE <- (rHo %*% solve(OMega1,
                         crossprod(aLpha.bar,
                                   crossprod(lm.fit(cbind(R0 %*% aLpha.ort,
                                                          R1),
                                                    R0)$res,
                                             lm.fit(cbind(R0 %*% aLpha.ort,
                                                          R1), R2)$res)
                         )
    ) +
      kAppa %*% solve(OMega2,
                      crossprod(aLpha.ort,
                                crossprod(R0,
                                          R1)
                      )
      )
    ) %*% Hmat
    tAu0 <- tAu
    tAu <- Hmat %*% t(matrix(solve(mB %x% mA + mD %x% mC,
                                   as.numeric(mE)
    ),
    ncol = ncol(Hmat)
    )
    )
    
    ll0 <- ll
    ll <- -0.5 * Time * c(determinant(crossprod(lm.fit(cbind(R2 %*% tAu %*% rHo +
                                                               R1 %*% pSi,
                                                             R1%*%tAu),
                                                       R0)$res) / Time
    )$modulus
    )
    ## cat(ll,'\t')
    ## if(ll < ll0 - 0.0001) stop('Likelihood is decreasing')
    if (abs((ll - ll0) / ll0) < 1e-12) break
  }
  aLpha.ort.bar <- aLpha.ort%*%solve(crossprod(aLpha.ort))
  oMega <- t(as.matrix(tmp1.r)[0:(p - r), , drop = FALSE])
  return(list(alpha = aLpha, tau = tAu, rho = rHo, psi = pSi, kappa = kAppa))
}

update.I1 <- function(obj, rank, ...) {
  # Update the specification of the I1 object
  # 
  #
  
  stopifnot(length(rank) %in% c(1, 2))
  if (length(rank) == 2 && rank[2] == 0) {
    obj$call[['rank']] <- rank[1]
    ans	<- eval(obj$call, envir = attr(obj$call, 'environment'))
  }
  else {
    obj[['call']][['rank']] <- rank
    ans <- eval(obj[['call']], envir = attr(obj$call, 'environment'))
  }
  return(ans)
}

VECM <- function(data, lags, dettrend, rank, season, exogenous, method = 'rrr') {
  # Function that creates a Vector Error Correction Model (VECM) object.
  #
  # Args:
  #   data: An xts object containing the data to be analyzed 
  #   lags: the number of lags in levels
  #   dettrend: specifies the deterministic trend structure
  #         - none
  #         - mean
  #         - drift
  #         - qtrend
  #   season: specified is whether inclusion of seasonal dummies is desired
  #   exogenous: an xts object with exogenous variables to be included 
  #       
  # Returns
  #   ft: an object of type I1
  #
  
  # Record the function call for later use
  cl <- match.call() # For the VECM call
  mf <- match.call() # For the CIModelFrame call
  
  m <- match(c('data', 'dummy', 'exogenous'), names(mf), 0L)
  if (any(missing(data), missing(lags), missing(dettrend))) stop('The arguments data, lags and dettrend are required')
  if (!inherits(data, 'xts')) stop('Data must be of class xts')
  if (any(is.na(data))) stop('Data contains NAs')
  if (!any(lags %in% seq(12))) stop ('The model must have lags between 1 and 12')
  if (!dettrend %in% c('none', 'mean', 'drift', 'qtrend', 'urmean')) stop('The determistic specification must be either none, mean or drift')
  
  # Check for inlcusion of exogenous variables as xts objects
  if (!missing(exogenous) && !is.null(exogenous)) {
    if (!inherits(exogenous, 'xts')) stop('Exogenous variables must be of class xts')
    if (sum(is.na(exogenous))) stop('Exogenous variables contains NAs') 
  }
  
  mf[[1L]] <- as.name("CIModelFrame")
  mf       <- eval(mf, parent.frame())
  
  if (missing(rank)) {
    ft      <- rrr(mf$Z0, mf$Z1, mf$Z2, ncol(mf$Z0))
    ft$lags <- lags
    ft      <- auxI1(ft)
  }
  else if (length(rank) == 1) {
    ft      <- rrr(mf$Z0, mf$Z1, mf$Z2, ncol(mf$Z0))
    ft$lags <- lags
    ft      <- setrankI1(ft, rank)
  }
  else if (length(rank) == 2) {
    ## if (!(cl[['dettrend']] %in% c('none', 'drift'))) stop('Deterministic specification must me either none og drift')
    if (lags < 2) stop('There must be at least tow lags in the model')
    mf$lags <- lags
    ft      <- setrankI2(mf, rank)
  }
  else stop('Illegal choice of rank')
  
  ft$call <- cl
  attr(ft$call, 'environment') <- parent.frame()
  ft$X    <- mf$X
  
  return(ft)
}

