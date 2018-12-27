#' @title Fit a specified cointegration model to data using reduced rank regression
#' 
#' @param data an xts object, time series data.
#' @param spec a ci_spec object, as returned by the function \code{ci_spec()} holding the
#'             specification information on the model
#' @param exogenous an xts object, time series data for included exogenous variables. These
#'                  could be more                   
#' @param H_alpha matrix, a matrix for imposing zero restrictions on alpha.
#' @param H_beta matrix, a design matrix with restrictions on the cointegration
#'        relations.
#' @param h_beta vector, the design vector with values for the restricted values
#'        in the cointegration relations
#' 
#' @export
#' @return a ci_fit object 

ci_fit <- function(data, 
                   spec = ci_spec(), 
                   H_alpha = NULL, 
                   H_beta = NULL,
                   h_beta = NULL,
                   exogenous = NULL) {
  
  if(is.null(H_alpha) && is.null(H_beta)) {
    dstruct <- ci_build_data_structures(data, spec, make_R = TRUE)
    estimates <- ci_estimate_unrestricted(dstruct$R0, dstruct$R1, spec$rank)
  }else {
    dstruct <- ci_build_data_structures(data, spec, make_S = TRUE)
    est_unres <- ci_estimate_unrestricted(dstruct$R0, dstruct$R1, spec$rank)
    estimates <- ci_switching_algorithm(est_unres, dstruct$S)
  }
  
  # Calculate the likelihood value ----
  
  # Put into a fit output vector ----
}











