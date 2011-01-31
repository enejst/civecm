#
# An R package for cointegration
# Author: Andreas Noack Jensen
# Date: April 2010
#
# Notes:
# The tracetest p-values has been stolen from McKinnon and I need to ask him
# in case I need to redistribute this code.
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
    if (ncol(ft$beta) > 0) {
        colnames(ft$beta) <- paste('beta(', 1:ncol(ft$beta), ')', sep = '')
        colnames(ft$alpha) <- paste('alpha(', 1:ncol(ft$alpha), ')', sep = '')
    }
    rownames(ft$beta) <- colnames(ft$Z1)
    rownames(ft$alpha) <- colnames(ft$Z0)
	
    if (ncol(ft$Z2)) {
        ft$Psi <- with(ft, t(lm.fit(Z2, Z0 - Z1 %*% tcrossprod(beta, alpha))$coef))
        colnames(ft$Psi) <- colnames(ft$Z2)
        ft$fitted.values <- with(ft, Z1 %*% tcrossprod(beta, alpha) +
                                 tcrossprod(Z2, Psi))
    }
    else ft$fitted.values <- with(ft, Z1 %*% tcrossprod(beta, alpha))
	
    ft$residuals <- ft$Z0 - ft$fitted.values
    colnames(ft$residuals) <- colnames(ft$Z0)
    class(ft$residuals) <- c('I1res', class(ft$residuals))

    class(ft) <- 'I1'
    return(ft)
}

auxI2 <- function(ft, r) {
    p0 <- ncol(ft$R0)
    p1 <- ncol(ft$R1)
    Time <- nrow(ft$R0)
    s1 <- if (is.null(ft$tau)) 0
    else ncol(ft$tau) - r
    s2 <- p0 - r - s1
    ft$beta <- ft$tau %*% ft$rho
    kSi <- -crossprod(ft$kappa, Null(ft$rho))
    eTa	<- crossprod(Null(ft$beta), ft$tau %*% Null(ft$rho))
    ft$delta <-	tcrossprod(Null(ft$tau)) %*% ft$psi
    tmpfit <- lm.fit(cbind(ft$R2 %*% ft$beta + ft$R1 %*% ft$delta, ft$R1 %*% ft$tau), ft$R0)
    tmpcoef <- t(matrix(tmpfit$coef, 2 * r + s1, p0))
    ft$alpha <-	tmpcoef[, seq_len(r), drop = FALSE]
    ft$zeta <- tmpcoef[, seq(r + 1, 2 * r + s1, len = r + s1), drop = FALSE]

    if (ncol(ft$Z3)) {
        ft$Psi <- t(lm.fit(ft$Z3, ft$Z0 - tcrossprod(ft$Z2 %*% ft$beta +
                                                     ft$Z1 %*% ft$delta, ft$alpha) -
                           ft$Z1 %*% tcrossprod(ft$tau, ft$zeta))$coef)
        colnames(ft$Psi) <- colnames(ft$Z3)
        ft$fitted.values <- with(ft, tcrossprod(Z2 %*% beta + Z1 %*% delta, alpha) +
                                 Z1 %*% tcrossprod(tau, zeta) +
                                 tcrossprod(Z3, Psi)
                                 )
    }
    else ft$fitted.values <- with(ft, tcrossprod(Z2 %*% beta + Z1 %*% delta, alpha) +
                                  ft$Z1 %*% tcrossprod(ft$tau, ft$zeta)
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

CIModelFrame	<- function(data, lags, dettrend, rank, season = FALSE, exogenous, dummy, ...) {
    stopifnot(inherits(data, 'matrix'))
    X <- data
    if (is.null(colnames(X))) colnames(X) <- paste('X', seq(ncol(X)), sep='.') 
    k <- lags
    Time <- nrow(X) - k
    if (dettrend == 0) {
        mTime <- matrix(rep(1, Time), , 1)
        colnames(mTime) <- 'constant'
    }
    else if(dettrend > 0) {
        mTime <- cbind(1, poly(seq_len(Time), dettrend, raw = T))
        colnames(mTime) <- c('constant', 'linear', 'quadratic', 'cubic', ifelse(dettrend > 3, NULL, paste('poly', seq(4, dettrend), sep = '.')))[seq_len(dettrend + 1)]
    }

    if (missing(exogenous)) W <- l.d.W <- l.d2.W <- NULL
    else {
        stopifnot(inherits(exogenous, 'matrix'))
        W <- exogenous
        stopifnot(nrow(W) == nrow(X))
        if (is.null(colnames(W))) colnames(W) <- paste('exo', seq_len(ncol(W)), sep = '.')
    }

    if (!season) Season <- NULL
    else {
        Season <- matrix(rep(c(1, rep(0, season - 2), -1), Time)[seq_len((season - 1) * Time)], Time, byrow = T)
        colnames(Season) <- paste('Season', seq_len(ncol(Season)), sep = '.')
    }

    if (missing(dummy)) Dummy <- NULL
    else {
        stopifnot(inherits(dummy, 'matrix'))
        stopifnot(nrow(dummy) == nrow(X))
        Dummy <- dummy[-seq_len(k),]
    }

    d.X	<- diff(X)
    colnames(d.X) <- paste('d', colnames(X), sep = '.')
    l.X <- X[seq_len(Time) + k - 1,]
    colnames(l.X) <- colnames(X)

    if (missing(rank) || length(rank) == 1) {
        l.d.X <- embed(d.X, k)
        colnames(l.d.X) <- c(paste('d', colnames(X), sep = '.'),
                                 if (k < 2) NULL
                                 else paste('l', rep(seq(k - 1), each = ncol(X)), '.d.', colnames(X), sep = ''))
    
        if (is.null(W)) l.d.W <- NULL
        else {
            l.d.W <- embed(diff(W), k)
            colnames(l.d.W) <- c(paste('d', colnames(W), sep = '.'),
                                 if (k < 2) NULL
                                 else paste('l', rep(seq(k - 1), each = ncol(W)), '.d.', colnames(W), sep = ''))
        }

        if (dettrend == -1) 
            Det1 <- Det2 <- NULL
        else if (dettrend == 0) {
            Det1 <- mTime[, 1, drop = F]
            Det2 <- NULL
        }
        else {
            Det1 <- mTime[, dettrend + 1, drop = F]
            Det2 <- mTime[, seq_len(dettrend), drop = F]
        }

        tmp <- list(X = X)
        tmp$Z0 <- l.d.X[, seq_len(ncol(X))]
        tmp$Z1 <- cbind(l.X, W[seq_len(Time) + k,], Det1)
        tmp$Z2 <- cbind(l.d.X[, -seq_len(ncol(X))], l.d.W, Season, Dummy, Det2)
        stopifnot(all(sapply(tmp, function(x) all(!is.na(x)))))
        class(tmp) <- 'model.frame.I1'
    }
    else if (length(rank) == 2) {
        if (k == 1) stop('In the I(2) model there must be at least two lags')

        l.d.X <- d.X[seq_len(Time) + k - 2,]
        colnames(l.d.X) <- paste('l', '.d.', colnames(X), sep = '')
        l.d2.X <- embed(diff(d.X), k - 1)
        colnames(l.d2.X) <-  c(paste('d2', colnames(X), sep = '.'),
                               if (k < 3) NULL
                               else paste('l', rep(seq(k - 2), each = ncol(X)), '.d2.', colnames(X), sep = ''))

        if (is.null(W)) d.W <- l.d2.W <- NULL
        else {
            d.W <- diff(W)[seq_len(Time) + k - 1,]
            colnames(d.W) <- paste('d', colnames(W), sep = '.')
            l.d2.W <- embed(diff(diff(W)), k - 1)
            colnames(l.d2.W) <- c(paste('d2', colnames(W), sep = '.'),
                                  if (k == 2) NULL
                                  else paste('l', rep(seq_len(k - 2), each = ncol(W)), '.d2.', colnames(W), sep =''))
        }
        
        if (dettrend == -1) 
            Det1 <- Det2 <- Det3 <- NULL
        else if (dettrend == 0) {
            Det2 <- mTime[, 1, drop = F]
            Det1 <- Det3 <- NULL
        }
        else if (dettrend == 1) {
            Det2 <- mTime[, 2, drop = F]
            Det1 <- mTime[, 1, drop = F]
            Det3 <- NULL
        }
        else {
            Det2 <- mTime[, dettrend + 1, drop = F]
            Det1 <- mTime[, dettrend, drop = F]
            Det3 <- mTime[, seq_len(dettrend - 1), drop = F]
        }

        tmp <- list(X = X)
        tmp$Z0 <- l.d2.X[,seq_len(ncol(X))]
        tmp$Z1 <- cbind(l.d.X, d.W, Det1)
        tmp$Z2 <- cbind(l.X, W[seq_len(Time) + k,], Det2)
        tmp$Z3 <- cbind(l.d2.X[,-seq_len(ncol(X))], l.d2.W, Season, Dummy, Det3)
        stopifnot(all(sapply(tmp, function(x) all(!is.na(x)))))
        class(tmp) <- 'model.frame.I2'
    }
	
    return(tmp)
}

## Fractional cointegration setup
## CIModelFrame	<-	function(data, lags, dettrend, model = c(d = 1, b = 1), poly.coint = FALSE, season = NULL, dummy = NULL, exogenous = NULL, ...)
# {
# 	X <- data
# 	if (is.null(colnames(X))) colnames(X) <- paste('X',seq(ncol(X)),sep='.')
# 	Time <- nrow(X) - lags
# 	d <- model[1]
# 	b <- model[2]

#     if (is.null(exogenous)) W <- l.d.W <- l.d2.W <- NULL
#     else if (any(model != 1)) stop('Model with exogenous varibles in fractional processes not handled yet.')
# 	else {
# 		W <- exogenous
# 		p.W <- ifelse(is.matrix(W), ncol(W), 1)
# 		l.d.W <- lag(diff(W), -(0:lags))
# 		if (is.null(colnames(W))) {
# 			dim(W) <- ifelse(is.null(dim(W)), c(Time, 1), dim(W))
# 			colnames(W) <- paste('exo', 1:p.W, sep = '.')
# 		}
# 		colnames(l.d.W) <- paste('.l.', rep(1:lags, each = p.W), 'd', colnames(W), sep ='')
# 		l.d2.W <- diff(l.d.W)
# 		colnames(l.d2.W) <- paste('l', rep(1:lags, each = p.W), '.d2.', colnames(W), sep ='')
# 	}

# 	if (!is.null(season)) {
# 		season <- ts(season * model.matrix(~factor(rep(1:season, length = Time)) - 1)
#                      - 1, end = end(X), freq = frequency(X))[, -season]
# 		colnames(season) <- paste('Csea', 1:ncol(season), sep = '.')
# 	}

# 	Dr	<- merge(W, switch(dettrend,
#                                   none = NULL,
#                                   mean = ts(c(rep(0, lags), rep(1, Time)), end = end(X),
#                                     freq = frequency(X)), drift = ts(c(rep(0, lags), seq(Time)),
#                                                             end = end(X), freq = frequency(X))))
# 	if(!is.null(dim(Dr))) {
# 		colnames(Dr) <- c(colnames(Dr), ifelse(dettrend != 'none', dettrend, NULL))
# 	}
# 	else if (!is.null(Dr)) {
# 		dim(Dr)	<- c(nrow(X), 1)
# 		colnames(Dr) <- dettrend
# 	}

# 	d.X <- fracdiff(X, d)
# 	colnames(d.X) <- paste('dd', colnames(X), sep = '.')
# 	bl.X <- X - fracdiff(X, b)
# 	if(b == d) db.bl.X	<- bl.X
# 	else if(b < d) db.bl.X	<- fracdiff(bl.X, d - b)
# 	else stop('Illega choice of b and d')
# 	colnames(db.bl.X) <- colnames(X)

# 	if (!poly.coint) {
# 		if(lags == 1) l.dd.bl.X	<- NULL
# 		else { 		
# 			l.dd.bl.X <- lag(fracdiff(bl.X, d), -(0:(lags - 2)))
# 			colnames(l.dd.bl.X) <- paste('l', rep(1:(lags - 1), each = ncol(X)), 
#                                          '.d.', colnames(X), sep = '')
# 		}
# 		Du <- merge(l.d.W, dummy, season
#                            , switch (dettrend,
#                                     none = NULL,
#                                     mean = NULL,
#                                     drift = ts(rep(1, Time), end = end(X), freq = frequency(X))))
# 		if (!is.null(dim(Du))) {
# 			colnames(Du) <- c(colnames(Du)[-length(colnames(Du))], ifelse(dettrend != 'none', 'const', NULL))
# 		}
# 		else if (!is.null(Du)) {
# 			dim(Du)	<- c(Time, 1)
# 			colnames(Du) <- 'const'
# 		}

# 		## Collect all Z's to get the timing and sample correct
# 		Ztmp <- na.omit(merge(d.X, db.bl.X, Dr, l.dd.bl.X, Du))
# 		colnames(Ztmp) <- unlist(sapply(list(d.X, db.bl.X, Dr, l.dd.bl.X, Du), colnames))

# 		tmp	<- list()
# 		tmp$X <- X
# 		tmp$Z0 <- Ztmp[,colnames(d.X), drop = F]
# 		tmp$Z1 <- Ztmp[,c(colnames(db.bl.X), colnames(Dr)), drop = F]
# 		tmp$Z2 <- if (!all(is.null(l.dd.bl.X), is.null(Du))) Ztmp[,c(colnames(l.dd.bl.X), 
#                                                                      colnames(Du)), drop = F]
#         else NULL
# 		class(tmp) <- 'model.frame.I1'
# 	}
# 	else {
# 		l.d.X <- lag(diff(X), -1)
# 		colnames(l.d.X) <- paste('l', '.d.', colnames(X), sep = '')
# 		d2.X <- diff(d.X)
# 		colnames(d2.X) <- paste('d2', colnames(X), sep = '.')
# 		if (lags == 1) stop('In the I(2) model there must be at least two lags')
# 		else if (lags == 2) l.d2.X <- NULL
# 		else {
# 			l.d2.X <- lag(diff(X, diff = 2), -(1:(lags - 2)))
# 			colnames(l.d2.X) <- paste('l', rep(1:(lags - 2), each = ncol(X)), '.d2.', colnames(X), sep = '')
# 		}
# 		if (!is.null(Dr)) {
# 			d.Dr <- diff(Dr)
# 			colnames(d.Dr) <- paste('d', colnames(Dr), sep = '.')			
# 		}
# 		else d.Dr <- NULL

# 		Du <- merge(l.d2.W, dummy, season)
# 		if (!is.null(dim(Du))) {
# 			colnames(Du) <- c(colnames(Du)[-length(colnames(Du))], ifelse(dettrend != 'none', 'const', NULL))
# 		}
# 		else if (!is.null(Du)) {
# 			dim(Du)	<- c(Time, 1)
# 			colnames(Du) <- 'trend'
# 		}

# 		## Collect all Z's to get the timing and sample correct
# 		Ztmp <- na.omit(merge(d2.X, l.d.X, d.Dr, bl.X, Dr, l.d2.X, Du))
# 		colnames(Ztmp) <- unlist(sapply(list(d2.X, l.d.X, d.Dr, bl.X, Dr, l.d2.X, Du), colnames))
	
# 		tmp	<- list()
# 		tmp$X <- X
# 		tmp$Z0 <- Ztmp[, colnames(d2.X), drop = F]
# 		tmp$Z1 <- Ztmp[, c(colnames(l.d.X), colnames(d.Dr)), drop = FALSE]
# 		tmp$Z2 <- Ztmp[, c(colnames(bl.X), colnames(Dr)), drop = FALSE]
# 		tmp$Z3 <- if (!all(is.null(l.d2.X), is.null(Du))) Ztmp[,c(colnames(l.d2.X), colnames(Du)), drop = FALSE]
#         else NULL
# 		class(tmp) <- 'model.frame.I2'
# 	}
	
# 	tmp$frequency <- frequency(Ztmp)
		
# 	return(tmp)
# }

## A function to print matrices of coefficient estimates with se og t.values
coefmatrix <- function(coef, stat1, digits = 3, sig.dig = c(TRUE, FALSE)) {
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

### Function borrowed from MASS but with here is the convention that the Null of a zero column matrix
### is the identity
Null <- function(M) {
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

    obj <- auxI1(obj)

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
    obj	<- auxI1(obj)
    
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
    obj	<- auxI2(obj, p0 - s1 - s2)
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
    else
        Z2.loc <- as.matrix(Z2)
    Z <- cbind(Z0.loc, Z1.loc, Z2.loc)
	
    Time <- nrow(Z)	
    p0 <- ncol(Z0.loc)
    p1 <- ncol(Z1.loc)
#	p2		<- ncol(Z2.loc)

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
    qrR1 <- qr(R1)
    tmpsvd <- svd(crossprod(qr.Q(qrR1), qr.Q(qr(R0))))
    ans	<- list()
    ans$values	<- tmpsvd$d^2
    ans$vectors	<- solve(qr.R(qrR1), tmpsvd$u * sqrt(Time))
    ans$R0 <- R0
    ans$R1 <- R1
    ans$S11 <- crossprod(R1) / Time
    ans$beta <- ans$vectors[, 1:r, drop = FALSE]
    ans$alpha <- t(lm.fit(ans$R1%*%ans$beta, ans$R0)$coef)
    ans$Z0 <- Z0
    ans$Z1 <- Z1
    ans$Z2 <- Z2
    ans$M <- M
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
    
    if (ncol(obj$Z3)) {
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
        initest <- rrr(diff(obj$X)[seq_len(Time) + obj$lags - 1,], obj$Z2, obj$Z1, r[1])
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
    ft <- auxI2(ft, r[1])
	
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
    if (ncol(object$Z2) && ncol(object$Z2) > (k - 1) * p)
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
    if (ncol(object$Z3) && ncol(object$Z3) > (k - 2) * p)
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
	
    object$t.alpha <- t.aLpha
    object$Pi <- PI
    object$t.Pi <- t.PI
    object$Psi <- object$Psi
    object$t.Psi <- t.PSi
    object$Omega <- OMega
    object$Sigma <- SIgma
    object$logLik <- ll
    object$rank <- rank
    object$p <- c(p0 = p0, p1 = p1, p2 = p2)
    
    class(object) <- 'summary.I1'
	
    return(object)
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

VECM <- function(data, lags, dettrend, rank, season, dummy, exogenous, method = 'rrr') {
    cl <- match.call()
    mf <- match.call()
    
    m <- match(c('data', 'dummy', 'exogenous'), names(mf), 0L)
    if (any(missing(data), missing(lags), missing(dettrend))) stop('The arguments data, lags and dettrend are required')
    if (any(is.na(data))) stop('Data contains NAs')
    if (!any(lags %in% seq(12))) stop ('The model must have lags between 1 and 12')

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

