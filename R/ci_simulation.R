#' @title Simulate data based on a ci_vecm model
#'
#' @param model a ci_vecm object
#' @param residuals xts object, residuals to generate data.
#'        primarily usefull in bootstrapping. Default is NULL
#'        in which case standard normal shocks are used.
#' 
#' @export
#' @return a matrix with with generated data

ci_simulate <- function(model, residuals = NULL) {
  
}