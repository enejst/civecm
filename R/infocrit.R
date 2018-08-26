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
