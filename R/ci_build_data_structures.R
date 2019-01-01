#' @title Build basic data structures for estimation and testing.
#' 
#' @description This function constructs the core data matrices used in
#'              in estimation and testing within the cointegrated vector
#'              error correction models. More specifically, the function
#'              builds the matrices Z0, Z1, Z2, R0, R1, S00, S10 and S11
#'              as these are described in 
#'              \insertCite{johansen_likelihood-based_1995}{civecm}.
#'
#' @param data an xts object, the series with the endogenous data.
#' @param spec a ci_spec object, the specification of the model.
#' @param make_R boolean, if TRUE; R matrices as discussed in Johansen(1996) 
#'        are constructed.
#' @param make_S boolean, if TRUE; S matrices as discussed in Johansen(1996)
#'        are constructed.
#' @param exogenous xts object, an time series with time index matching that
#'        of data. 
#'
#' @importFrom Rdpack reprompt
#'
#' @references{
#'   \insertRef{johansen_likelihood-based_1995}{civecm}
#' }
#'
#' @export
#' @return a list with lists of Z, R and S matrices.

ci_build_data_structures <- function(data, 
                                     spec = ci_spec(), 
                                     make_R = TRUE, 
                                     make_S = TRUE, 
                                     exogenous = NULL) {
  
  Z <- list(Z0 = NULL, Z1 = NULL, Z2 = NULL)
  R <- list(R0 = NULL, R1 = NULL)
  S <- list(S00 = NULL, S01 = NULL, S11 = NULL)
  
  p <- ncol(data)
  q <- spec$lags
  r <- spec$rank
  n <- nrow(data)
  m <- n - q
  det <- spec$det_spec
  dates <- zoo::index(data)
  pe <- ncol(exogenous)
  
  Z$Z1 <- xts::lag.xts(data, 1)[-(1:q),]
  Z$Z0 <- data[-(1:q),] - Z$Z1
  if(q > 1) {
    Z$Z2 <- xts::lag.xts(data - xts::lag.xts(data, 1), seq_len(q-1))[-(1:q),]
  }
  
  Z <- ci_add_deterministics_to_Z(det, q, Z)
  
  if(!is.null(exogenous)) {
    Z <- ci_add_exogenous_to_Z(q, Z, exogenous)
  }
    
  if(make_R || make_S) {
    if(q == 1) {
      R$R0 <- Z$Z0
      R$R1 <- Z$Z1
    }else {
      R$R0 <- stats::lm.fit(Z$Z2, Z$Z0)
      R$R1 <- stats::lm.fit(Z$Z2, Z$Z1)
    }
    
    if(make_S) {
      S$S00 <- (1/n) * crossprod(R$R0, R$R0)
      S$S01 <- (1/n) * crossprod(R$R0, R$R1)
      S$S11 <- (1/n) * crossprod(R$R1, R$R1)
    }
  }
  data_structures <- list(Z = Z, R = R, S = S)
}