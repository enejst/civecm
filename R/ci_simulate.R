#' @title Simulate data based on a ci_vecm model
#'
#' @param model a ci_vecm object, the model to simulate from.
#' @param nobs integer, the number of simulated observations. Ignored if
#'        shocks or exogenous are provided.
#' @param initial xts object, a q x p xts object with initial values for the 
#'        endogenous variables. If not provided (default is NULL), 
#'        a matrix of zeroes q x p zeroes is used. 
#' @param shocks xts object, nobs x p shocks to generate data.
#'        Primarily usefull in bootstrapping. Default is NULL
#'        in which case standard normal shocks are used together with 
#'        Omega from the model object.
#' @param exogenous xts object, a (nobs + q) x pe matrix of exogenous variables 
#'        that can be considered fixed. Here, q denotes the number of lags in 
#'        the model. Mostly relevant for bootstrapping in a model with 
#'        exogenous regressors. If shocks are provided, the number of rows 
#'        in this matrix must match the number of rows in shocks plus 
#'        the number of model lags.
#' 
#' @export
#' @return an xts object with with generated data

ci_simulate <- function(model,
                        nobs = 100,
                        initial = NULL,
                        shocks = NULL, 
                        exogenous = NULL) {

  # Get model data ----
  p <- nrow(model$alpha)
  lags <- model$spec$lags
  detspec <- model$spec$deterministic
  alpha <- model$alpha
  beta <- model$beta
  tGamma <- ifelse(is.null(model$Gamma), NULL, t(model$Gamma))
  tPhi <- ifelse(is.null(model$Phi), NULL, t(model$Phi))
  
  # Input validation ---- 
  # We check that if shocks and exogenous are both provided, 
  # the number of observations in both match
  if(!is.null(shocks)) {
    if(!is.null(exogenous)) {
      stopifnot((nrows(shocks) + lags) == nrows(exogenous))
    }
  }
  
  # Verify the initial values input
  stopifnot((class(initial) == "xts") || is.null(initial))
  
  # Construct holders for output and simulation data ----
  if(class(exogenous) == "xts") {
    time_idx <- zoo::index(exogenous)[-(1:lags)]
    
    if(is.null(initial)) {
      initial <- matrix(0, lags, p)
      initial <- xts::xts(initial, zoo::index(exogenous)[(1:lags)])
    }
  }else if(class(shocks == "xts")) {
    time_idx <- zoo::index(shocks)
  }else {
    time_idx <- seq(Sys.Date() - nobs, Sys.Date(), by = "days")
  }
  
  if(!is.null(shocks)) {
    nobs <- nrow(shocks) + lags
  }else {
    Omega <- model$Omega
    shocks <- xts::xts(randn(nobs, p) %*% Omega, time_idx)
  }
  
  if(is.null(initial)) {
    per <- periodicity(time_idx)$scale
    initial <- matrix(0, lags, p)
    initial <- xts::xts(initial, 
                        seq(time_idx[1]-(lags + 1), time_idx - 1, by = per))
  }
  X <- xts::xts(matrix(NA, nobs, p), time_idx)
  X <- rbind(initial, X)
  
  D <- switch(
    detspec,
    "ci_constant" = xts::xts(matrix(1, nobs, 1), time_idx),
    "constant" = xts::xts(matrix(1, nobs, 1), time_idx),
    "ci_trend" = xts::xts(cbind(matrix(1, nobs, 1), matrix(1:nobs, nobs, 1)), time_idx),
    "trend" = xts::xts(cbind(matrix(1, nobs, 1), matrix(1:nobs, nobs, 1)), time_idx),
    "none" = NULL,
    stop("Unknown deterministic specification"))
  
  if(!is.null(exogenous)) {
    pe <- ncol(E)
    E <- exogenous
    LE <- xts::lag.xts(E, 1)
    DE <- E - LE
    DEE <- xts::lag.xts(DE, seq_len(lags - 1))
  }
  
  # Run the simulation loop ---- 
  for(i in 1 + lags:(nobs + lags)) {
    X[i,] <- X[i - 1,] + 
             X[i - 1] %*% tcrossprod(beta[1:p,], alpha) + shocks[i - lags,]
    
    if(lags > 1) {
      DX_1 <- diff(X[(i - lags):(i - 1)], 1)
      DXX <- tail(xts::lag.xts(DX_1, seq_len(lags - 2)), 1)
      X[i,] <- X[i,] + DXX %*% tGamma
    }
    
    X[i,] <- X[i,] + switch(
      detspec,
      "ci_constant" = D %*% tcrossprod(beta[p,], alpha),
      "constant" = D %*% tcrossprod(beta[p,], alpha) + tPhi[1,],
      "ci_trend" = D %*% tcrossprod(beta[p + 0:1,], alpha) + tPhi[1,],
      "trend" = D %*% (tcrossprod(beta[p + 0:1,], alpha) + tPhi[1:2,]),
      "none" = 0,
      stop("Unknown deterministic specification")
    )
    
    if(!is.null(exogenous)) {
      X[i,] <- X[i,] + LE[i - 1,] %*% tcrossprod(beta[p + 1:pe,], alpha)
      if(lags > 1) {
        iDEE <- tail(xts::lag.xts(DEE, seq_len(lags - 1)), 1)
        itPhi <- switch(
          detspec, 
          "ci_constant" = tPhi, 
          "constant" = tPhi[2:pe,], 
          "ci_trend" = tPhi[2:pe,],
          "trend" = tPhi[2:pe,],
          "none" = NULL,
          stop("Unknown deterministic specification")
        )
        X[i,] <- X[i,] + iDEE * itPhi
      }
    }
  }
  X <- X[-(1:lags),]
}