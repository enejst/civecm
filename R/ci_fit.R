#' @title Fit a specified cointegration model to data using reduced rank regression
#' 
#' @param spec a ci_spec object, as returned by the function \code{ci_spec()} holding the
#'             specification information on the model
#' @param data an xts object, time series data.
#' @param exogenous an xts object, time series data for included exogenous variables. These
#'                  could be more 
#' @param H_alpha list, a list of H matrices for implosing linear restrictions on alpha
#'                as described in Juselius(2006, ch. 9). 
#' @param H_beta list, a list of H matrices for imposing linear restrictions on beta as
#'               in Juselius(2006, ch. 8)
#' 
#'
#' @export
#' @return a ci_fit object 

ci_fit <- function(spec = ci_spec(), 
                   data, 
                   H_alpha = NULL, 
                   H_beta = NULL,
                   exogenous = NULL) {
  
  # Specifications ----
  p <- ncol(data)
  q <- spec$lags
  r <- spec$rank
  n <- nrow(data)
  m <- n - q
  det <- spec$det_spec
  dates <- zoo::index(data)
  pe <- ncol(exogenous)
  
  # Build the data structures ----
  LX <- lag(data, 1)[-(1:q),]
  DX <- (data - lag(data,1))
  
  # Construct the lagged differences ----
  LDX <- lag(DX, 1:(q-1))[-(1:(q-1)),]
  Z0 <- LDX[,1:p]
  LDX <- LDX[,-(1:p)]
    
  # Add deterministics to the cointegration relations ----
  if(det == "none") {
    Z1 <- LX
    Z2 <- LDX
  }else if(det == "ci_constant") {
    constant <- xts::xts(1, dates)
    Z1 <- merge(Z1, constant)
    Z2 <- LDX
  }else if(det == "constant") {
    constant <- xts::xts(1, dates)
    Z1 <- merge(Z1, constant)
    Z2 <- merge(LDX, constant)
  }else if(det == "ci_trend") {
    constant <- xts::xts(1, dates)
    trend <- xts::xts(seq_len(n-q), dates)
    Z1 <- merge(Z1, constant)
    Z1 <- merge(Z1, trend)
    Z2 <- merge(LDX, constant)
  }else if(det == "trend") {
    constant <- xts::xts(1, dates)
    trend <- xts::xts(seq_len(n-q), dates)
    Z1 <- merge(Z1, constant)
    Z1 <- merge(Z1, trend)
    Z2 <- merge(LDX, constant)
    Z2 <- merge(Z2, trend)
  }
  
  # Construct basic H_alpha and H_beta is not provided
  
  
}











