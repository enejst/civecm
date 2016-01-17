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
# Date 09/12/2015
# Development taken over by Emil Nejstgaard
#
#

anova.I1 <- function(object, ..., df) {
    if(length(objects <- list(object, ...)) == 1) stop('Nothing to test against')
    test <- 2 * (logLik(objects[[1]]) - logLik(objects[[2]]))
    return(data.frame(Value = test,
                      Test = paste('ChiSq(', df, ')', sep = ''),
                      p.value = 1 - pchisq(test, df)
                      )
           )
}

auxI1 <- function(ft) {
    if (ncol(ft$beta) > 0)
    {
      colnames(ft$beta)  <- paste('beta(', 1:ncol(ft$beta), ')', sep = '')
      colnames(ft$alpha) <- paste('alpha(', 1:ncol(ft$alpha), ')', sep = '')
    }
    rownames(ft$beta)  <- colnames(ft$Z1)
    rownames(ft$alpha) <- colnames(ft$Z0)
    
    if (!is.null(ft$Z2))
    {
      tmpfit           <- with(ft, lm.fit(Z2, Z0 - Z1 %*% tcrossprod(beta, alpha)))
      ft$Psi           <- t(tmpfit$coef)
      colnames(ft$Psi) <- colnames(ft$Z2)
      ft$fitted.values <- xts(tmpfit$fitted.values, index(ft$Z0))
      ft$residuals     <- xts(tmpfit$residuals, index(ft$Z0))
    }
    else
    {
      ft$fitted.values <- with(ft, as.xts(Z1 %*% tcrossprod(beta, alpha), index(Z0)))
      ft$residuals     <- with(ft, as.xts(Z0 - fitted.values, index(Z0)))
    }
  
    ft$Omega               <- crossprod(ft$residuals) / nrow(ft$residuals)
    colnames(ft$residuals) <- colnames(ft$Z0)
    class(ft$residuals)    <- c('I1res', class(ft$residuals))
  
    class(ft) <- 'I1'
    return(ft)
}

auxI2 <- function(ft, r) {
    p0   <- ncol(ft$R0)
    p1   <- ncol(ft$R1)
    Time <- nrow(ft$R0)
    s1   <- if (is.null(ft$tau)) 0
    else ncol(ft$tau) - r
    s2       <- p0 - r - s1
    ft$beta  <- ft$tau %*% ft$rho
    kSi      <- -crossprod(ft$kappa, Null(ft$rho))
    eTa	     <- crossprod(Null(ft$beta), ft$tau %*% Null(ft$rho))
    ft$delta <-	tcrossprod(Null(ft$tau)) %*% ft$psi
    tmpfit   <- lm.fit(cbind(ft$R2 %*% ft$beta + ft$R1 %*% ft$delta, ft$R1 %*% ft$tau), ft$R0)
    tmpcoef  <- t(matrix(tmpfit$coef, 2 * r + s1, p0))
    ft$alpha <-	tmpcoef[, seq_len(r), drop = FALSE]
    ft$zeta  <- tmpcoef[, seq(r + 1, 2 * r + s1, len = r + s1), drop = FALSE]

    if (!is.null(ncol(ft$Z3))) {
        ft$Psi <- t(lm.fit(ft$Z3, ft$Z0 - tcrossprod(ft$Z2 %*% ft$beta +
                                                     ft$Z1 %*% ft$delta, ft$alpha) -
                           ft$Z1 %*% tcrossprod(ft$tau, ft$zeta))$coef)
        colnames(ft$Psi) <- colnames(ft$Z3)
        ft$fitted.values <- with(ft, xts(tcrossprod(Z2 %*% beta + Z1 %*% delta, alpha) +
                                         Z1 %*% tcrossprod(tau, zeta) +
                                         tcrossprod(Z3, Psi),
                                         index(Z0)
                                         )
                                 )
    }
    else ft$fitted.values <- with(ft, xts(tcrossprod(Z2 %*% beta + Z1 %*% delta, alpha) +
                                          ft$Z1 %*% tcrossprod(ft$tau, ft$zeta),
                                          index(Z0)
                                          )
                                  )

    if (s1 != 0) {
        ft$beta1 <- Null(ft$beta) %*% eTa
        ft$alpha1 <- Null(ft$alpha) %*% kSi
        rownames(ft$beta1) <- colnames(ft$R2)
        colnames(ft$beta1) <- paste('beta1(', 1:ncol(ft$beta1), ')', sep = '')
        rownames(ft$alpha1)	<- colnames(ft$R0)
        colnames(ft$alpha1) <- paste('alpha1(', 1:ncol(ft$alpha1), ')', sep = '')
    }
    if (s2 != 0) {
        ft$beta2 <- Null(ft$beta) %*% Null(eTa)
        ft$alpha2 <- Null(ft$alpha) %*% Null(kSi)
        rownames(ft$beta2) <- colnames(ft$R2)
        colnames(ft$beta2) <- paste('beta2(', 1:ncol(ft$beta2), ')', sep = '')
        rownames(ft$alpha2) <- colnames(ft$R0)
        colnames(ft$alpha2) <- paste('alpha2(', 1:ncol(ft$alpha2), ')', sep = '')
    }

    if (r > 0) {
        colnames(ft$beta) <- paste('beta(', 1:ncol(ft$beta), ')', sep = '')
        colnames(ft$alpha) <- paste('alpha(', 1:ncol(ft$alpha), ')', sep = '')
        rownames(ft$beta) <- c(colnames(ft$R2))
    }

    ft$residuals <- ft$Z0 - ft$fitted.values
    colnames(ft$residuals) <- colnames(ft$Z0)
    class(ft$residuals) <- c('I2res', class(ft$residuals))

    class(ft) <- c('I2', 'I1')

    return(ft)
}

betaNorm <- function(obj, ...) UseMethod('betaNorm')

betaNorm.default <- function(obj, norm.matrix, ...) {
    if (is.logical(norm.matrix)) obj <- t(t(obj) / obj[norm.matrix])
    else obj <- obj %*% solve(norm.matrix)
    return(obj)
}

betaNorm.I1 <- function(obj, norm.matrix, ...) {
    obj$beta <- if (missing(norm.matrix))
        betaNorm(obj$beta, qr.R(qr(obj$R1%*%obj$beta, LAPACK=T))/sqrt(nrow(obj$R0)))
    else betaNorm(obj$beta, norm.matrix)
    obj$alpha <- t(lm.fit(obj$R1%*%obj$beta, obj$R0)$coef)
    return(obj)
}

boot <- function(bootfct, obj, reps, wild = TRUE, cluster = NULL) {
    reslen <- nrow(obj$residuals)
    sf <- function(dummy) {
        newcall	<- obj$call
        if (wild)

            newcall[['data']] <- simulate(obj, res = obj$residuals * rnorm(reslen))
        else
            newcall[['data']] <- simulate(obj,
                                          res = obj$residuals[sample(seq(reslen),
                                          replace = TRUE),]
                                          )
        return(bootfct(eval(newcall, envir = attr(obj$call, 'environment'))))
    }
    if (is.null(cluster)) {
        if (suppressWarnings(require('pbapply'))) ans <- pblapply(seq(reps), sf)
        else ans <- lapply(seq(reps), sf)
    }
    else {
        clusterEvalQ(cluster, library(civecm))
        if (wild) clusterSetupRNG(cluster, type = 'RNGstream')
        ans <- parLapply(cluster, seq(reps), sf)
    }
    return(ans)
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
  #       - urmean
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
  if (missing(exogenous)) W <- l.W <- l.d.W <- l.d2.W <- NULL
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
    Z2     <- merge(l.d.X, l.d.W, l.d.detTrend, Season)
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
    # A function to print matrices of coefficient estimates with se og t.values
  
    m.coef <- as.matrix(coef)
    m.stat1 <- as.matrix(stat1)
    mchar <- max(nchar(cbind(m.coef, m.stat1)))
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
    maxi <- apply(m.out, 2, nchar)
    li <- apply(m.out, 2, function(x) regexpr('.', x, fixed = TRUE))
    ri <- maxi - li
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
  #       - urmean
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
  else if(dettrend == "drift") {
    # Construct xts variable with a constant and polynomial trends
    detTrend               <- xts(cbind(poly(seq_len(length(dateIndex)), 1, raw = T)), dateIndex)
    l.detTrend             <- lag(detTrend)
    colnames(l.detTrend)   <- c('ci-trend')
    l.d.detTrend           <- diff(lag(detTrend))
    colnames(l.d.detTrend) <- "constant"
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
  } else if(dettrend == "urmean") {
    detTrend             <- xts(rep(1, length(dateIndex)), dateIndex)
    l.detTrend           <- NULL
    l.d.detTrend         <- lag(detTrend)[-seq_len(lags)]
    colnames(l.d.detTrend) <- 'constant'
    l.d2.detTrend        <- NULL
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

### Calculates the eigenvalues of the companion matrix of the VAR. The argument
### 'roots' can be either a vector of ranks for which the eigenvalues should be
### computed, the default value 'max' which will give the eigenvalues for all
### possible cointegration ranks. The values will be calculated without any
### restrictions on the parameters. Finally it is possible parse the value
### 'specific' and the function will give the eigenvalues for the specifik model
### with restrictions impoced.

eigen <- function(object, ...) UseMethod('eigen')

eigen.matrix <- function(object, ...) base:::eigen(object, ...)

eigen.I1 <- function(object, ...) {
    out <- eigen(companion(VAR(object)), ...)
    class(out)	<- 'eigenCompanion'
    return(out)
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
  p2     <- ifelse(is.missing(obj$Z2) || is.null(obj$Z2), NULL, ncol(obj$Z2))
  r      <- obj$rank
  G      <- ifelse(is.missing(obj$Ha), diag(p0*r), obj$Ha)   
  H      <- ifelse(is.missing(obj$Hb), diag(p1*r), obj$Hb)
  alpha  <- obj$alpha
  beta   <- obj$beta
  Omega  <- obj$Omega
  iOmega <- solve(Omega)
  T      <- nrow(obj$Z1)
  S11    <- obj$M[(p0+1):(p0+p1),(p0+1):(p0+p1)]
  
  hessian_11 <- t(G) %*% kronecker(iOmega, crossprod(beta,S11) %*% beta) %*% G
  hessian_12 <- t(G) %*% kronecker(iOmega %*% alpha, crossprod(beta,S11)) %*% H 
  hessian_22 <- t(H) %*% kronecker(crossprod(alpha,iOmega) %*% S11) %*% H
  
  hessian <- T * rbind(cbind(hessian_11,hessian_12),
                       cbind(hessian_21,hessian_22))
  
  urHessian_11 <- kronecker(iOmega, crossprod(beta,S11) %*% beta)
  urHessian_12 <- kronecker(iOmega %*% alpha, crossprod(beta,S11))
  urHessian_22 <- kronecker(crossprod(alpha,iOmega) %*% S11)
  
  urHessian <- T * rbind(cbind(urHessian_11,urHessian_12),
                         cbind(urHessian_21,urHessian_22))
  
  sdPsi   <- sqrt(diag(hessian[ncol(G),ncol(G)]))
  sdPhi   <- sqrt(diag(hessian[ncol(H),ncol(H)]))
  sdAlpha <- t(matrix(sdPsi,nrow=r,ncol=p0))
  sdBeta  <- matrix(sdPhi,nrow=p0,ncol=r)
  
  covarList           <- list(hessian=hessian)
  covarList$iHessian  <- solve(hessian)
  covarList$urHessian <- urHessian
  covarList$sdPI      <- sqrt(diag(urHessian))
  covarList$sePi      <- covarList$sdPI / sqrt(T)
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

    obj	<- aux.I1(obj)
    class(obj) <- 'I1gls'

    return(obj)
}

infocrit <- function(obj) {
    if (!inherits(obj, 'I1')) stop('Object must be of typee either I1 or summary.I1')

    sm <-summary(obj)
    Time <- nrow(sm$residuals)
    logOm <- as.double(determinant(sm$Omega)$modulus)
    sc <- as.double(logOm + (sm$p['p2'] * sm$p['p0'] +
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

residuals.I1	<- function(object, ...) {
    tmp	<- object$residuals
    return(tmp)
}

restrictAlpha <- function(obj, A.matrix, type = 'restr') {
    if(!inherits(A.matrix, 'matrix')) stop('Restriction matrix must be a matrix')
    else if (type == 'restr') {
        H.bar <- A.matrix %*% solve(crossprod(A.matrix))
        H.ort <- Null(A.matrix)
        tmp.rrr <- rrr(obj$R0 %*% H.bar, obj$R1, obj$R0 %*% H.ort, ncol(obj$beta))
        obj$H.alpha <- A.matrix
        obj$alpha <- A.matrix %*% tmp.rrr$alpha
        obj$beta <- tmp.rrr$beta
    }
    else if (type == 'known') {
        a.bar <- A.matrix %*% solve(crossprod(A.matrix))
        a.ort <- Null(A.matrix)
        a.bar.ort <- Null(a.bar)
        if (ncol(A.matrix) == ncol(obj$alpha)) {
            obj$alpha <- A.matrix
            obj$beta <- solve(crossprod(obj$R1), crossprod(obj$R1, obj$R0 %*% a.bar))
        }
        else {
            tmp.rrr <- rrr(obj$R0 %*% a.bar.ort,
                           obj$R1,
                           ,
                           ncol(obj$beta) - ncol(A.matrix)
                           )
            tmp.X <- cbind(obj$R1,
                           obj$R0 %*%
                           a.bar.ort - obj$R1 %*%
                           tcrossprod(tmp.rrr$beta, tmp.rrr$alpha)
                           )
            tmp.fit <- solve(crossprod(tmp.X), crossprod(tmp.X, obj$R0 %*% a.bar))
            obj$alpha <- cbind(A.matrix, a.ort %*% tmp.rrr$alpha)
            obj$beta <- cbind(tmp.fit[1:ncol(obj$Z1), ], tmp.rrr$beta)
        }
    }

    obj <- aux.I1(obj)

    return(obj)
}

restrictBeta <- function(obj, H.matrix) {
    if (length(H.matrix) == 1) H.matrix	<- H.matrix[[1]]
    if (inherits(H.matrix, 'matrix')) {
        temp.rrr <- rrr(obj$Z0, obj$Z1 %*% H.matrix, obj$Z2, ncol(obj$beta))
        obj$H.beta <- H.matrix
        obj$beta <- H.matrix %*% temp.rrr$beta
        ## obj$beta <- betaNorm(obj$beta, obj$S11)
        obj$alpha <- t(lm.fit(obj$R1%*%obj$beta, obj$R0)$coef)
    }
    else if (!inherits(H.matrix, 'list'))
        stop('H.matrix must be either a matrix or a list of matrices.')
    else if (length(H.matrix) != ncol(obj$beta))
        stop('The number of elements in H.matrix must match the cointegration rank')
    else {
        for (i in seq_len(ncol(obj$beta))) {
            obj$beta[, i] <- obj$beta %*%
                geigen(crossprod(obj$beta, H.matrix[[i]]) %*%
                       solve(crossprod(H.matrix[[i]]),
                             crossprod(H.matrix[[i]], obj$beta)
                             ),
                       crossprod(obj$beta)
                       )$vectors[, 1]
        }

        i <- 1
	loglik0 <- logLik(obj)
        iter <- 1

        while(1) {
            temp.rrr <- rrr(obj$R0,
                            obj$R1 %*% H.matrix[[i]],
                            obj$R1 %*% obj$beta[, -i, drop = F],
                            1
                            )

            obj$beta[, i] <- H.matrix[[i]] %*% temp.rrr$beta

            if (i == ncol(obj$beta)) {
                ##obj$beta		<- betaNorm(obj$beta, obj$S11)
                tmp.lm <- lm.fit(obj$R1%*%obj$beta, obj$R0)
                obj$alpha <- t(tmp.lm$coef)
                obj$residuals <- tmp.lm$residuals
                loglik1 <- logLik(obj)
                if (loglik1<loglik0 && iter != 1) {
                    print(c(loglik0, loglik1))
                    stop('Loglik is decreasing. Try to impose identifying restrictions')
                }
                if (abs((loglik0 - loglik1) / loglik0) < 1e-20) {
                    cat('Convergence in', iter, 'iterations.\n')
                    break
                }
                loglik0	<- loglik1
                i <- 1
                iter <- iter + 1
            }
            else i <- i + 1
        }
    }
    obj	<- aux.I1(obj)

    return(obj)
}

restrictLongRun <- function(obj,Hb,hb,Ha,alpha_i,beta_i) {
  # Function that restricts both subelements of the long run matrix
  # The function applies the general Boswijk/Doornik principle for
  # imposing linear restrictions on the long run parameters.
  # See Boswijk and Doornik (2004) : "Identifying, Estimating and
  # Testing restricted cointegration systems"
  #
  # Args
  #   obj: an I1 object with restricted parameter estimates
  #   Hb: matrix with linear restrictions on vec(beta)
  #   hb: a vector of corresponding fixed parameter values
  #   Ha: matrix with linear restrictions on vec(alpha)
  #   alpha_i: a matrix with initial values for alpha
  #   beta_i: a matrix with initial values for beta
  #
  # Returns
  #   obj: an object of type I1  
  #
  
  p0            <- ncol(obj$Z0)
  p1            <- ncol(obj$Z1)
  p2            <- ifelse(is.null(obj$Z2), 0, ncol(obj$Z2))
  p             <- p0
  r             <- ncol(alpha_i)
  p1            <- nrow(beta_i)
  obj_i         <- obj
  obj_i$alpha   <- alpha_i
  obj_i$beta    <- beta_i
  obj_i         <- auxI1(obj_i)
  obj_i$Hb      <- Hb
  obj_i$hb      <- hb
  obj_i$Ha      <- Ha
  S11           <- obj_i$M[(p0+1):(p0+p1),(p0+1):(p0+p1),drop=F]
  S01           <- obj_i$M[(1:p0),(p0+1):(p0+p1),drop=F]
  S00           <- obj_i$M[(1:p0),(1:p0),drop=F]
  alpha_prime_i <- t(alpha_i)  
  PI_hat_LS     <- solve(S11,t(S01))
  Omega_i       <- (S00 - S01 %*% tcrossprod(beta_i,alpha_i) 
                    - t(S01 %*% tcrossprod(beta_i,alpha_i))
                    + tcrossprod(alpha_i,beta_i) %*% S11 %*% tcrossprod(beta_i,alpha_i))
  invOmega_i    <- solve(Omega_i)
  
  loglik_i      <- logLik(obj_i)
  loglik_1      <- loglik_i - 1
  
  i <- 1
  while( (loglik_i - loglik_1) > 1e-6 && i < 10000) {
    # Reset loglikelihood
    loglik_1 <- loglik_i
    
    # Estimate alpha
    psi_i          <- solve(crossprod(Ha, kronecker(invOmega_i,t(beta_i) %*% S11 %*% beta_i) %*% Ha), 
                            crossprod(Ha, kronecker(invOmega_i,t(beta_i) %*% S11))) %*% vec(PI_hat_LS)
    valpha_prime_i <- Ha %*% psi_i
    alpha_prime_i  <- matrix(valpha_prime_i,ncol=p)
    alpha_i        <- t(alpha_prime_i)
    
    # Estimate beta
    phi_i   <- solve(crossprod(Hb, kronecker(crossprod(alpha_i,invOmega_i %*% alpha_i),S11) %*% Hb),
                     crossprod(Hb, kronecker(crossprod(alpha_i,invOmega_i),S11))) %*% (vec(PI_hat_LS) - kronecker(alpha_i,diag(p1)) %*% hb)
    vbeta_i <- Hb %*% phi_i + hb
    beta_i  <- matrix(vbeta_i,ncol=r)
    
    # Estimate Omega
    Omega_i    <- (S00 - S01 %*% tcrossprod(beta_i,alpha_i) - t(S01 %*% tcrossprod(beta_i,alpha_i))
                   + tcrossprod(alpha_i,beta_i) %*% S11 %*% tcrossprod(beta_i,alpha_i))
    invOmega_i <- solve(Omega_i)
    
    # Calculate likelihood
    obj_i$alpha <- alpha_i
    obj_i$beta  <- beta_i
    
    obj_i    <- auxI1(obj_i)
    loglik_i <- logLik(obj_i)
    i        <- i + 1
  }
  obj <- auxI1(obj_i);
  return(obj)
}

restrictTau <- function(obj, Hmat) {
    if (!inherits(obj, 'I2')) stop('Object must have class I2')
    if (nrow(Hmat) != ncol(obj$R1))
        stop('Restriction matrix has wrong number of rows')
    if (!(ncol(Hmat) %in% seq(ncol(obj$R1))))
        stop('Restriction matrix has wrong number of columns')
    Time <- nrow(obj$R0)
    p0 <- ncol(obj$R0)
    p1 <- ncol(obj$R1)
    s1 <- ncol(obj$beta1)
    s2 <- ncol(obj$beta)

    if (p0 - s1 - s2 > 0)
        ans <- .Fortran('tauswitch',
                        as.double(obj$R0),
                        as.double(obj$R1),
                        as.double(obj$R2),
                        as.double(Hmat),
                        time = as.integer(Time),
                        as.integer(p0),
                        as.integer(p1),
                        as.integer(p0 - s1 - s2),
                        as.integer(s2),
                        as.integer(ncol(Hmat)),
                        tau = obj$tau,
                        rho = matrix(0.0, s2 + s1, p0 - s1 - s2),
                        psi = matrix(0.0, p1, s2),
                        kappa = matrix(0.0, s2 + s1, p0 - s2),
                        as.integer(10000), 1e-9, integer(1))
    else
        ans <- restrictBeta(update(obj, s1), Hmat)

    obj$tau <- ans$tau
    obj$rho <- ans$rho
    obj$psi <- ans$psi
    obj$kappa <- ans$kappa
    obj	<- aux.I2(obj, p0 - s1 - s2)
    return(obj)
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
  #    p0 <- ncol(Z0.loc)
  #    p1 <- ncol(Z1.loc)
  #	   p2	<- ncol(Z2.loc)
  
  M <- crossprod(Z) / Time
  
  if (is.null(Z2.loc)) {
    R0 <- Z0
    R1 <- Z1
  }
  else {
    R0 <- lm.fit(Z2, Z0)$residuals
    R1 <- lm.fit(Z2, Z1)$residuals
  }
  ##	ans <- .Fortran('rrr', as.double(R0), as.double(R1), as.integer(Time),
  ## as.integer(p0), as.integer(p1), values = double(p1),
  ## vectors = matrix(0.0, p1, p1))
  
  ## The following is to be implemented in Fortran
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

    obj <- aux.I1(obj)
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

simulate.I1 <- function(object, nsim, seed, res, init, ...) {
    stopifnot(!missing(res))

    p <- nrow(object$alpha)
    r <- ncol(object$alpha)
    k <- object$lags

    A <- VAR(object)

    if(missing(init)) init <- object$X[seq_len(k), ]

    eps <- as.matrix(res) +
        tcrossprod(object$Z1[, seq_len(ncol(object$Z1) - p) + p, drop = F] %*%
                   object$beta[seq_len(ncol(object$Z1) - p) + p, , drop = F],
                   object$alpha
                   )
    if (!is.null(ncol(object$Z2)) && ncol(object$Z2) > (k - 1) * p)
        eps <-  eps + tcrossprod(object$Z2[, seq_len(ncol(object$Z2) - (k - 1) * p) +
                                           (k - 1) * p, drop = F],
                                 object$Psi[, seq_len(ncol(object$Z2) - (k - 1) * p) +
                                            (k - 1) * p, drop = F]
                                 )
    attributes(eps) <- attributes(res)
    ans <- rbind(init, rmfilter(eps, A, init))

    return(ans)
}

simulate.I2 <- function(object, nsim, seed, res, init, ...) {
    stopifnot(!missing(res))

    p <- nrow(object$alpha)
    r <- ncol(object$alpha)
    k <- object$lags

    A <- VAR(object)

   	if(missing(init)) init <- object$X[seq_len(k), ]

    eps <- as.matrix(res) +
        tcrossprod(object$Z2[, seq_len(ncol(object$Z2) - p) + p, drop = F] %*%
                   object$beta[seq_len(ncol(object$Z2) - p) + p, , drop = F] +
                   object$Z1[, seq_len(ncol(object$Z1) - p) + p, drop = F] %*%
                   object$delta[seq_len(ncol(object$Z1) - p) + p, , drop = F],
                   object$alpha
                   ) + object$Z1[, seq_len(ncol(object$Z1) - p) + p, drop = F] %*%
                       tcrossprod(object$tau[seq_len(ncol(object$Z1) - p) + p,
                                             , drop = F],
                                  object$zeta
                                  )
    if (!is.null(ncol(object$Z3)) && ncol(object$Z3) > (k - 2) * p)
      eps <- eps + tcrossprod(object$Z3[, seq_len(ncol(object$Z3) - (k - 2) * p) +
                                        (k - 2) * p, drop = F],
                              object$Psi[, seq_len(ncol(object$Z3) - (k - 2) * p) +
                                         (k - 2) * p, drop = F]
                              )
    attributes(eps) <- attributes(res)

    ans <- rbind(init, rmfilter(eps, A, init))

    return(ans)
}

summary.I1 <- function(object, ...) {
    p0 <- ncol(object$Z0)
    p1 <- ncol(object$Z1)
    p2 <- ifelse(is.null(object$Z2), 0, ncol(object$Z2))
    Time <- nrow(object$Z0)
    rank <- ncol(object$beta)

    PI <- if (ncol(object$beta) != 0) tcrossprod(object$alpha, object$beta)

    OMega <- crossprod(object$residuals) / Time
    SIgma <- if (!is.null(object$Z2)) {
        if (!is.null(PI)) object$M[-(1:(p0 + p1)), -(1:(p0 + p1)), drop = F] -
            object$M[-(1:(p0 + p1)), (p0 + 1):(p0 + p1), drop = F] %*%
                object$beta %*%
                    ginv(crossprod(object$beta,
                                   object$M[(p0 + 1):(p0 + p1), (p0 + 1):(p0 + p1)]) %*%
                         object$beta
                         ) %*%
                             crossprod(object$beta,
                                       object$M[(p0 + 1):(p0 + p1), -(1:(p0 + p1))]
                                       )
    }
    else NULL

    if (!is.null(PI)) {
        se.aLpha <- sqrt(diag(OMega) %o%
                         diag(ginv(crossprod(object$beta,
                                             object$S11) %*% object$beta)) / Time
                         )
		t.aLpha	<- object$alpha / se.aLpha

        se.PI <- sqrt(diag(OMega) %o%
                      diag(object$beta %*%
                           tcrossprod(ginv(crossprod(object$beta, object$S11) %*%
                                           object$beta), object$beta)) / Time
                      )
        t.PI <- PI / se.PI
    }
    else se.aLpha <- t.aLpha <- se.PI <- t.PI <- NULL

    ll <- -nrow(object$Z0) * 0.5 * c(determinant(OMega)$modulus) -
        Time * p0 * 0.5 * (log(2 * pi) + 1)

    if (!is.null(object$Z2) & !is.null(PI)) {
        se.PSi <- sqrt(diag(OMega) %o% diag(solve(SIgma)) / Time)
        t.PSi <- object$Psi / se.PSi
    }
    else se.PSi <- t.PSi <- NULL
    colnames(PI)

    tmp	<- object['call']
    tmp$beta <- object$beta
    tmp$alpha <- object$alpha
    tmp$t.alpha <- t.aLpha
    tmp$Pi <- PI
    tmp$t.Pi <- t.PI
    tmp$Psi <- object$Psi
    tmp$t.Psi <- t.PSi
    tmp$Omega <- OMega
    tmp$Sigma <- SIgma
    tmp$eigenvalues <-object$values
    tmp$logLik <- ll
    tmp$rank <- rank
    tmp$p <- c(p0 = p0, p1 = p1, p2 = p2)
    tmp$lags <- object$lags
    tmp$residuals <- object$residuals

    class(tmp) <- 'summary.I1'

    return(tmp)
}

summary.I2 <- function(object, ...) {
    p0 <- ncol(object$Z0)
    p1 <- ncol(object$Z1)
    p2 <- ncol(object$Z2)
    p3 <- ifelse(is.null(object$Z3), 0, ncol(object$Z3))
    Time <- nrow(object$Z0)
    object$Omega <- crossprod(object$residuals) / Time
    object$logLik <- -nrow(object$Z0) * 0.5 * c(determinant(object$Omega)$modulus) -
        Time * p0 * 0.5 * log(2 * pi * exp(1))

    class(object) <- 'summary.I2'

    return(object)
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
            H <- diag(nrow(obj$alpha))[, -i]
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

tracetest <- function(fit, type, a.p.value, bootstrap, reps, cluster) UseMethod('tracetest')

tracetest.I1 <- function(fit, type = 'I1', a.p.value = TRUE, bootstrap = FALSE, reps = 499, cluster = NULL) {
    if (type == 'I1') {
        p0 <- ncol(fit$Z0)
        ## Johansen 1995, second edition
        ## MacKinnon et al 1998. JAE. UDEN TILLADELSE...
        ans <- data.frame(p.r = p0:1)
        ans$r <- 0:(p0 - 1)
        ans$eig.value <- fit$values[1:p0]
        ans$Trace <- -nrow(fit$Z0) * rev(cumsum(rev(log(1 - fit$values[1:p0]))))
        if (a.p.value) {
            if (fit$call[['dettrend']] == -1) dettrend <- 'none'
            else if (fit$call[['dettrend']] == 0) dettrend <- 'mean'
            else if (fit$call[['dettrend']] == 1) dettrend <- 'drift'
            else stop('Deterministic specification not covered')
            ans	<- cbind(ans,
                         tracemckinnon(ans$p.r,
                                       'trace',
                                       dettrend,
                                       c(0.05, 0.01)
                                       ),
                         tracemckinnon(ans$p.r,
                                       'trace',
                                       dettrend,
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

update.I1 <- function(object, rank, ...) {
    stopifnot(length(rank) %in% c(1, 2))
    if (length(rank) == 2 && rank[2] == 0) {
        object$call[['rank']] <- rank[1]
        ans	<- eval(object$call, envir = attr(object$call, 'environment'))
    }
    else {
        object[['call']][['rank']] <- rank
        ans <- eval(object[['call']], envir = attr(object$call, 'environment'))
    }
    return(ans)
}

VAR <- function(object, ...) UseMethod('VAR')

VAR.I1 <- function(object, ...) {
    p <- ncol(object$X)
    A <- array(0, c(p, p, object$lags))
    A[,,1] <- with(object, diag(p) + tcrossprod(alpha, beta[seq_len(p),, drop = FALSE]))
    for (i in seq_len(object$lags - 1)) {
        A[, , i] <- A[, , i] + object$Psi[, seq_len(p) + (i - 1) * p]
        A[, , i + 1] <- A[, , i + 1] - object$Psi[, seq_len(p) + (i - 1) * p]
    }

    class(A) <- 'VAR'
    return(A)
}

VAR.I2 <- function(object, ...) {
    p <- ncol(object$X)

    A <- array(0, c(p, p, object$lags))

    A[, , 1] <- with(object, 2 * diag(p) +
                     tcrossprod(alpha, tau[seq_len(p), ] %*% rho + delta[seq_len(p),]) +
                     tcrossprod(zeta, tau[seq_len(p), ]))
    A[, , 2] <- with(object, -(diag(p) + tcrossprod(alpha, delta[seq_len(p), ]) +
                  tcrossprod(zeta, tau[seq_len(p), ])))

    for (i in seq_len(max(object$lags - 2, 0))) {
        A[, , i] <- A[, , i] + object$Psi[, seq_len(p) + (i - 1) * p]
        A[, , i + 1] <- A[, , i + 1] - 2 * object$Psi[, seq_len(p) + (i - 1) * p]
        A[, , i + 2] <- A[, , i + 2] + object$Psi[, seq_len(p) + (i - 1) * p]
    }

    class(A) <- 'VAR'
    return(A)
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
  #         - urmean
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
  if (!dettrend %in% c('none', 'mean', 'drift','qtrend','urmean')) stop('The determistic specification must be either none, mean or drift')
  
  # Check for inlcusion of exogenous variables as xts objects
  if (!missing(exogenous)) {
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


