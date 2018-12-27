#' @title The estimation switching algorithm to estimate the model under
#'        general linear restrictions
#'        
#' @description This function applies the method outlined in e.g. Boswijk(1995).
#'        
#' @param unres_fit list, a list with unrestricted estimates for alpha and beta.
#' @param S list, a set of product moment matrices, S00, S01, S11 as 
#'        described in e.g. Boswijk(1995) or Johansen(1995).
#' @param H_alpha matrix, a matrix for imposing zero restrictions on alpha.
#' @param H_beta matrix, a design matrix with restrictions on the cointegration
#'        relations.
#' @param h_beta vector, the design vector with values for the restricted values
#'        in the cointegration relations
#' 
#' @export
#' @return a list with the matrices alpha_hat, beta_hat and omega_hat

ci_switching_algorithm <- function(unres_fit, 
                                   S, 
                                   H_alpha = NULL, 
                                   H_beta = NULL, 
                                   h_beta = NULL) {
  
  # Dimensions ----
  p <- nrow(unres_fit$alpha_hat)
  r <- ncol(unres_fit$alpha_hat)
  px <- nrow(unres_fit$beta_hat)
  
  # Check for H_alpha and construct standard if not found ----
  if(is.null(H_alpha)) {
    H_alpha <- diag(p * r)
  }
  
  # Check for H_beta ----
  if(is.null(H_beta)) {
    H_beta <- diag(px * r)
    h_beta <- matrix(0, px * r, 1)
  }
  
  # Get initial values from unrestricted model ----
  v_beta_hat_ur <- as.vector(unres_fit$beta_hat)
  i_v_phi <- solve(tcrossprod(H_beta,H_beta),
                   H_beta %*% as.vector(unres_fit$beta_hat))
  i_v_gam <- solve(tcrossprod(H_alpha,H_alpha),
                   H_alpha %*% as.vector(unres_fit$alpha_hat))
  
  # Do switching algorithm ----
  
}