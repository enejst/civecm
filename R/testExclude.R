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
