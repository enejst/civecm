#' @title Calculate the residuals in the model
#'
#' @param alpha matrix, a p x r matrix with error correction paramters
#' @param beta matrix, a d x r matrix with long run relation, where d
#'        is at least p and further depends on the deterministic 
#'        specification and exogenous variables
#' @param R0 matrix, a T x p matrix with concentrated residuals
#' @param R1 matrix, a T x d matrix with concentrated residuals
#' 
#' @return a T x p matrix with model residuals

ci_calculate_residuals <- function(alpha, beta, R0, R1) {
  
  e <- R0 - R1 %*% tcrossprod(beta, alpha)
  
}