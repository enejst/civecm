
testResiduals <- function(obj) {
  sm <- summary(obj)
  ans <- list()
  
  ans$cor <- sm$Omega / sqrt(diag(sm$Omega) %o% diag(sm$Omega))
  ans$se <- sqrt(diag(sm$Omega))
  
  ans$info <- infocrit(obj)
  ans$autocor <- testAutocor(obj)
  ans$normal <- testNormal(residuals(obj))
  ans$arch <- testARCH(obj)
  
  return(ans)
}
