### Calculates the eigenvalues of the companion matrix of the VAR. The argument
### 'roots' can be either a vector of ranks for which the eigenvalues should be
### computed, the default value 'max' which will give the eigenvalues for all
### possible cointegration ranks. The values will be calculated without any
### restrictions on the parameters. Finally it is possible parse the value
### 'specific' and the function will give the eigenvalues for the specifik model
### with restrictions impoced.

eigen <- function(object, ...) UseMethod('eigen')

eigen.matrix <- function(object, ...) base:::eigen(object, ...)

eigen.I1 <- function(object, ...) {
  out <- eigen(companion(VAR(object)), ...)
  class(out)	<- 'eigenCompanion'
  return(out)
}