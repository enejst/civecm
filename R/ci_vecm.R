#' @title Create a cointegrated Vector Error Correction model object. 
#'
#' @param spec a ci_spec object, a ci_spec object defining the model specification.
#' @param alpha a matrix, the p x r matrix of error correction parameters.
#' @param beta a matrix, the matrix of long run relations which is at least p x r.
#' @param Omega a matrix, the positive semidefinite covariance matrix of the residuals.
#' @param Gamma a matrix, the short run parameters.
#' @param Phi a matrix, the coefficients for the exogenous variables outside the 
#'        cointegration relations.
#'
#' @export 
#' @return ci_model object

ci_vecm <- function(spec, alpha, beta, Omega, Gamma = NULL, Phi = NULL) {
    
  model <- list(
    spec = spec,
    alpha = alpha,
    beta = beta,
    Gamma = Gamma,
    Omega = Omega,
    Phi = Phi)
  
  class(model) <- "ci_vecm"
  
  model
}