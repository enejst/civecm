betaNorm <- function(obj, ...) UseMethod('betaNorm')

betaNorm.default <- function(obj, norm.matrix, ...) {
  if (is.logical(norm.matrix)) obj <- t(t(obj) / obj[norm.matrix])
  else obj <- obj %*% solve(norm.matrix)
  return(obj)
}

betaNorm.I1 <- function(obj, norm.matrix, ...) {
  if (missing(norm.matrix)) {
    obj$beta <- betaNorm(obj$beta, qr.R(qr(obj$R1%*%obj$beta, LAPACK=T))/sqrt(nrow(obj$R0)))
  }
  else {
    obj$beta <- betaNorm(obj$beta, norm.matrix)
  }
  obj$alpha <- t(lm.fit(obj$R1%*%obj$beta, obj$R0)$coef)
  return(obj)
}