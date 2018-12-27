#' @title Calculate the error covariance matrix Omega
#'
#' @param residuals, a T x p matrix of model residuals
#'
#' @return a matrix with the estimated error covariance

ci_calculate_omega <- function(residuals) {
  Omega <- crossprod(residuals) / nrow(residuals)
}