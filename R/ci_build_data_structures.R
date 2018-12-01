ci_build_data_structures <- function(spec, data, make_M = TRUE, make_R = TRUE, make_S = TRUE, exogenous = NULL) {
  
  
  # Output List ----
  Z <- list(Z0 = NULL, Z1 = NULL, Z2 = NULL)
  M <- list(M02 = NULL, M12 = NULL, M22 = NULL)
  R <- list(R0 = NULL, R1 = NULL)
  S <- list(S00 = NULL, S01 = NULL, S11 = NULL)
  
  # Build data structures ----
  if(make_Z = TRUE || make_M = TRUE || make_R == TRUE || make_S = TRUE) {
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
    if(q > 1) Z$Z2 <- lag(Z0, 1:(q-1))[-(1:(q-1)),]
    
    Z <- ci_add_deterministics_to_Z(det, q, Z)
    
    if(!is.null(exogenous)) {
      Z <- ci_add_exogenous_to_Z(det, q, Z, exogenous)
    }
    
    if(make_M == TRUE || make_R == TRUE || make_S == TRUE) {
      if(q > 1) {
        if(M == TRUE) {
          M$M02 <- (1/n) * crossprod(Z0, Z2)
          M$M12 <- (1/n) * crossprod(Z1, Z2)
          M$M22 <- (1/n) * crossprod(Z2, Z2)
        }
      }
      # Construct R matrices (residuals) ----
      if(make_R == TRUE || make_S == TRUE) {
        if(q == 1) {
          R$R0 <- Z0
          R$R1 <- Z1
        }else {
          M22iZ2 <- solve(M22, Z2)
          R$R0 <- Z$Z0 - M02 %*% M22iZ2
          R$R1 <- Z$Z1 - M12 %*% M22iZ2
        }
      
        if(make_S == TRUE) {
          S$S00 <- (1/n) * crossprod(R0, R0)
          S$S01 <- (1/n) * crossprod(R0, R1)
          S$S11 <- (1/n) * crossprod(R1, R1)
        }
      }
    }
  }
  
  data_structures <- list(Z = Z, M = NULL, R = NULL, S = NULL)
  
  
}