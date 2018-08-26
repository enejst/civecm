auxI1 <- function(ft)
{
  if (ncol(ft$beta) > 0)
  {
    colnames(ft$beta) <- paste('beta(', 1:ncol(ft$beta), ')', sep = '')
    colnames(ft$alpha) <- paste('alpha(', 1:ncol(ft$alpha), ')', sep = '')
  }
  rownames(ft$beta) <- colnames(ft$Z1)
  rownames(ft$alpha) <- colnames(ft$Z0)
  
  if (!is.null(ft$Z2))
  {
    tmpfit <- with(ft, lm.fit(Z2, Z0 - Z1 %*% tcrossprod(beta, alpha)))
    ft$Psi <- t(tmpfit$coef)
    colnames(ft$Psi) <- colnames(ft$Z2)
    ft$fitted.values <- xts(tmpfit$fitted.values, index(ft$Z0))
    ft$residuals <- xts(tmpfit$residuals, index(ft$Z0))
  }
  else
  {
    ft$fitted.values <- with(ft, as.xts(Z1 %*% tcrossprod(beta, alpha), index(Z0)))
    ft$residuals <- with(ft, as.xts(Z0 - fitted.values, index(Z0)))
  }
  
  colnames(ft$residuals) <- colnames(ft$Z0)
  class(ft$residuals) <- c('I1res', class(ft$residuals))
  
  class(ft) <- 'I1'
  return(ft)
}

