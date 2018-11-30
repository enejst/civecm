ci_add_deterministics_to_Z <- function(det_spec, q, Z) {
  
  Z0 <- Z$Z0
  Z1 <- Z$Z1
  Z2 <- Z$Z2
  
  if(det == "none") {
    Z1 <- LX
  }else if(det == "ci_constant") {
    constant <- xts::xts(1, dates)
    Z1 <- merge(Z1, constant)
  }else if(det == "constant") {
    constant <- xts::xts(1, dates)
    Z1 <- merge(Z1, constant)
    Z2 <- merge(Z2, constant)
  }else if(det == "ci_trend") {
    constant <- xts::xts(1, dates)
    trend <- xts::xts(seq_len(n-q), dates)
    .  <- merge(Z1, constant)
    Z1 <- merge(., trend)
    Z2 <- merge(Z2, constant)
  }else if(det == "trend") {
    constant <- xts::xts(1, dates)
    trend <- xts::xts(seq_len(n-q), dates)
    .  <- merge(Z1, constant)
    Z1 <- merge(., trend)
    if(q > 1) {
      .  <- merge(Z2, constant)
      Z2 <- merge(., trend)
    }
  
  list(Z0 = Z0, Z1 = Z1, Z2 = Z2)
}