#' @title Create a specification for a cointegrated vector autoregression
#' 
#' @description Builds a specification for a standard cointegrated VAR model.
#'              It includes specifications of lag structure, cointegration rank
#'              and basic deterministic specifications.
#'              
#' @param lags numeric integer, number of lags in levels
#' @param rank numeric integer, the cointegration rank
#' @param det_spec character, the base deterministic specification. The ones
#'        supported are the five given in Juselius(2006).
#'        \itemize{
#'          \item none - no deterministics
#'          \item ci_constant - constant restricted to the cointegration relations
#'          \item constant - unrestricted constant
#'          \item ci_trend - unrestricted constant and trend restricted to the
#'                cointegration relations
#'          \item trend - unrestricted constant and trend
#'        }
#' @param seasonals boolean, if TRUE seasonal dummies corresponding to the
#'        data frequency will be included, monthly and quarterly are supported.
#'
#' @export
#' @return creates an object of type ci_spec

ci_spec <- function(lags = 2, 
                    rank = 1,
                    det_spec = "constant", 
                    seasonals = FALSE) {
  
  spec <- list(
    lags = lags,
    rank = rank,
    det_spec = det_spec,
    seasonals = seasonals
  )
  class(spec) <- "ci_spec"
  
  spec
}