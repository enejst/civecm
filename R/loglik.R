logLik.I1 <- function(object, const = TRUE, ...) {
  OMega <- with(object, crossprod(residuals) / nrow(Z0))
  ll <- -nrow(object$Z0) * 0.5 * c(determinant(OMega)$modulus)
  if (const) ll <- ll - prod(dim(object$Z0)) * 0.5 * (log(2 * pi) + 1)
  names(ll) <- 'logLik'
  return(ll)
}