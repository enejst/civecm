#' @title Fit a specified cointegration model to data using reduced rank regression
#' 
#' @param spec a ci_spec object, as returned by the function \code{ci_spec()} holding the
#'             specification information on the model
#' @param data an xts object, time series data.
#' @param exogenous an xts object, time series data for included exogenous variables. These
#'                  could be more                   
#' @param alpha_res list, a list containing a vector h_alpha and a matrix H_alpha giving
#'        the general linear restrictions as in Boswijk(1995) as well as initial values
#'        for the freely varying parameter in alpha, i_alpha. 
#' @param beta_res list, a list containing a vector h_beta and a matrix H_beta giving
#'        the design for general linear restrictions as in Boswijk(1995) as well as
#'        initial values for the freely varying parameters in beta called, i_beta.
#' 
#' @export
#' @return a ci_fit object 

ci_fit <- function(spec = ci_spec(), 
                   data, 
                   alpha_res = list(H_alpha = NULL, h_alpha = NULL, i_alpha = NULL), 
                   beta_res = list(H_beta = NULL, h_beta = NULL, i_beta = NULL),
                   exogenous = NULL) {
  
  # Construct data structures
  data_structures <- ci_build_data_structures(spec, data, Z = TRUE)
  
  # Construct basic H_alpha and H_beta if not provided ----
  if(is.null(alpha_res$H_alpha) && is.null(beta_res$H_beta)) {
    estimates <- ci_estimate_unrestricted()
  }else {
    estimates <- ci_switching_algo()
  }
  
  # Calculate the likelihood value ----
  
  # Put into a fit output vector ----
}











