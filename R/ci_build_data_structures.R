#' @title 
#'
#' @param data an xts object, the series with the endogenous data.
#' @param spec a ci_spec object, the specification of the model.
#' @param make_R boolean, if TRUE R matrices as discussed in Johansen(1996) 
#'        are constructed.
#' @param make_S boolean, if TRUE S matrices as discussed in Johansen(1996)
#'        are constructed.
#'        
#' @export
#' @return a list with Z, R and S matrices, only Z matrices are madatory.

ci_build_data_structures <- function(data, 
                                     spec = ci_spec(), 
                                     make_R = TRUE, 
                                     make_S = TRUE, 
                                     exogenous = NULL) {
  
  
  # Output List ----
  Z <- list(Z0 = NULL, Z1 = NULL, Z2 = NULL)
  R <- list(R0 = NULL, R1 = NULL)
  S <- list(S00 = NULL, S01 = NULL, S11 = NULL)
  
  # Build data structures ----
  if(make_Z = TRUE || make_R == TRUE || make_S = TRUE) {
    p <- ncol(data)
    q <- spec$lags
    r <- spec$rank
    n <- nrow(data)
    m <- n - q
    det <- spec$det_spec
    dates <- zoo::index(data)
    pe <- ncol(exogenous)
    
    Z$Z0 <- (data - lag(data,1))
    Z$Z1 <- lag(data, 1)[-(1:q),]
    if(q = 1) Z$Z2 = NULL
    if(q > 1) Z$Z2 <- xts::lag.xts(Z0, 1:(q-1))[-(1:(q-1)),]
    
    Z <- ci_add_deterministics_to_Z(det, q, Z)
    
    if(!is.null(exogenous)) {
      Z <- ci_add_exogenous_to_Z(det, q, Z, exogenous)
    }
    
    # Construct R matrices (residuals) ----
    if(make_R == TRUE || make_S == TRUE) {
      if(q == 1) {
        R$R0 <- Z0
        R$R1 <- Z1
      }else {
        R$R0 <- lm.fit(Z2, Z0)
        R$R1 <- lm.fit(Z2, Z1)
      }
      
      # Construct S matrices ----
      if(make_S == TRUE) {
        S$S00 <- (1/n) * crossprod(R0, R0)
        S$S01 <- (1/n) * crossprod(R0, R1)
        S$S11 <- (1/n) * crossprod(R1, R1)
      }
    }
  }
  
  data_structures <- list(Z = Z, R = R, S = S)
}