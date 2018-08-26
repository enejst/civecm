rstandard.I1 <- function(model, ...) {
  return(model$residuals / matrix(sqrt(diag(summary(model)$Omega)),
                                  nrow(model$residuals),
                                  ncol(model$residuals),
                                  byrow = T
  )
  )
}
