#' @title Calculates the eigenvalues of the companion matrix of the VAR.
#' @description  The values will be calculated without any restrictions on the parameters. 
#'               
#' @param roots atomic vector or character, a vector of ranks or the default value 'max' 
#'              which will give the eigenvalues for all ossible cointegration ranks.
#'              It is possible parse the value 'specific' and the function will give the 
#'              eigenvalues for the specific model with restrictions impoced.
#' 
#' @export
#' @return an object of class eigenCompanion

eigen <- function(object, ...) UseMethod('eigen')

eigen.matrix <- function(object, ...) base:::eigen(object, ...)

eigen.I1 <- function(object, ...) {
  out <- eigen(companion(VAR(object)), ...)
  class(out)	<- 'eigenCompanion'
  return(out)
}
