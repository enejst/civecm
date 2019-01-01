#' @title Add exogenous variables to Z1 and Z2
#'
#' @param lags, the number of lags in the model 
#' @param Z list, the list with the Z matrices
#' @param exogenous xts object with the exogenous variables
#'
#' @export
#' @return the list with Z0, Z1 and Z2

ci_add_exogenous_to_Z <- function(lags, Z, exogenous) {
  
  Z1 <- Z$Z1
  Z2 <- Z$Z2
  
  E1 <- xts::lag.xts(exogenous)
  Z1 <- merge(Z1, E1)
  Z$Z1 <- Z1
  
  if(lags > 1) {
    DE <- exogenous - E1
    DEE <- xts::lag.xts(DE, seq_len(lags - 1))[-(1:(lags-1)), ]
    Z2 <- merge(Z2, DEE)
    Z$Z2 <- Z2
  }
  
  Z
}