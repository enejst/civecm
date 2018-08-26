testUnit <- function(obj) {
  if(!inherits(obj, 'I1')) stop('Object must have I1')
  
  r.max <- ncol(obj$X) - 1
  
  values <- p.vals <- matrix(NA, r.max, ncol(obj$X))
  tmp <- matrix('', 2 * r.max, nrow(obj$alpha))
  
  for (r in seq_len(r.max)) {
    fb <- update(obj, r)
    for (i in seq_len(nrow(obj$alpha))) {
      H <- diag(nrow(obj$alpha))[, i, drop = FALSE]
      tmp <- anova.I1(fb,
                      restrictAlpha(fb, H, ifelse(r != 1, 'known', 'restr')),
                      df = ncol(obj$X) - r
      )
      values[r, i] <- tmp[1, 1]
      p.vals[r, i] <- tmp[1, 3]
    }
  }
  
  colnames(values) <- colnames(obj$X)
  out <- list(df = ncol(obj$X) - 1:r.max, value = values, p.value = p.vals)
  class(out) <- 'autotest.I1'
  return(out)
}
