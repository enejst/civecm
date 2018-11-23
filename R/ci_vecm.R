#' @title Create a cointegrated Vector Error Correction model object. 
#'
#' @param spec a ci_spec object, a ci_spec object defining the model specification
#' @param alpha a matrix, the p x r matrix of error correction parameters.
#' @param beta a matrix, the matrix of long run relations which is at least p x r.
#' @param Gamma a matrix, the short run parameters
#' @param Omega a matrix, the positive semidefinite covariance matrix of the residuals.
#' @param Lambda a matrix, the short run parameters for potential exogenous variables.
#'
#' @export 
#' @return ci_model object

ci_vecm <- function(spec, alpha, beta, Gamma, Omega, Lambda = NULL) {
  
  
  
}