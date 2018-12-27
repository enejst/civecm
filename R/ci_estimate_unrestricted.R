#' @title Estimate the unrestricted model by Reduced Rank Regression
#'
#' @param R0 matrix, as returned by the function ci_build_data_structures
#' @param R1 matrix, as returned by the function ci_build_data_structures
#' @param r integer, the number of cointegration relations
#'
#' @return a list with an estimate for alpha and beta.

ci_estimate_unrestricted <- function(R0, R1, r) {
  svdR0 <- svd(R0)
  svdR1 <- svd(R1)
  svd01 <- svd(crossprod(svdR0$u, svdR1$u))
  values	<- svd01$d^2
  vectors	<- svd01$v * sqrt(Time)
  for (i in 1:length(svdR1$d))
  {
    if (svdR1$d[i] > 2e-16)
      ans$vectors[i,] <- vectors[i,]/svdR1$d[i]
    else
      ans$vectors[i,] <- 0
  }
  vectors <- svdR1$v%*%ans$vectors
  beta_hat <- vectors[, 1:r, drop = FALSE]
  alpha_hat <- t(lm.fit(R1%*%beta_hat, R0)$coef)
  
  results <- list(alpha = alpha_hat, beta = beta_hat)
}