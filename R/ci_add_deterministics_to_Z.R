#' @title Add deterministic terms to Z
#' 
#' @param det_spec character, the deterministic specification
#' @param q numeric, the number of lags in the model
#' @param Z list, a list with Z0, Z1, Z2
#'
#' @return list with Z matrices updated with deterministics

ci_add_deterministics_to_Z <- function(det_spec, q, Z) {
  
  Z0 <- Z$Z0
  Z1 <- Z$Z1
  Z2 <- Z$Z2
  
  constant <- xts::xts(1, dates)
  trend <- xts::xts(seq_len(n-q), dates)
  
  Z1 <- switch(
    det,
    "none" = LX,
    "ci_constant" = merge(Z1, constant),
    "constant" = merge(Z1, trend),
    "ci_trend" = merge(merge(Z1, constant), trend),
    "trend" = merge(merge(Z1, constant), trend),
    stop("Unknown deterministic specification.")
  )
  
  Z2 <- switch(
    det,
    "none" = NULL,
    "ci_constant" = NULL,
    "constant" = merge(Z2, constant),
    "ci_trend" = merge(Z2, constant),
    "trend" = { 
      if(q > 1) {
        .  <- merge(Z2, constant)
        Z2 <- merge(., trend)
      }else {
        Z2 <- merge()
      }
      Z2
    },
    stop("Unknown deterministic specification")
  )
  
  list(Z0 = Z0, Z1 = Z1, Z2 = Z2)
}