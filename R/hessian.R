hessian <- function(obj) {
  # Calculates the covariance matrix for the parameter
  # estimates given restrictions imposed through the function
  # restrictLongRun
  #
  # Args:
  #   obj: an object of type I1
  #
  # Return
  #   covarList: a list with covariance matrices and 
  #
  
  p0     <- ncol(obj$Z0)
  p1     <- ncol(obj$Z1)
  p2     <- ifelse(is.null(obj$Z2), NULL, ncol(obj$Z2))
  G      <- if(is.null(obj$Ha)) { diag(p0*r) } else {obj$Ha}    
  H      <- if(is.null(obj$Hb)) { diag(p1*r) } else {obj$Hb}
  alpha  <- obj$alpha
  beta   <- obj$beta
  Omega  <- obj$Omega
  iOmega <- solve(Omega)
  r      <- ncol(alpha)
  T      <- nrow(obj$Z1)
  S11    <- crossprod(obj$R1,obj$R1) / T
  
  hessian_11 <- t(G) %*% kronecker(iOmega, crossprod(beta,S11) %*% beta) %*% G
  hessian_12 <- t(G) %*% kronecker(iOmega %*% alpha, crossprod(beta,S11)) %*% H 
  hessian_22 <- t(H) %*% kronecker(crossprod(alpha,iOmega) %*% alpha, S11) %*% H
  
  hessian <- T * rbind(cbind(hessian_11,hessian_12),
                       cbind(t(hessian_12),hessian_22))
  iHessian <- solve(hessian)
  
  sdPsi   <- if(ncol(G)==1) {
    sqrt(iHessian[1:ncol(G),1:ncol(G)])
  }
  else {
    sqrt(diag(iHessian[1:ncol(G),1:ncol(G)]))}
  sdPhi   <- if(ncol(H)==1) {
    sqrt(iHessian[(ncol(G)+1):(ncol(G)+ncol(H)),(ncol(G)+1):(ncol(G)+ncol(H))])} 
  else {
    sqrt(diag(iHessian[(ncol(G)+1):(ncol(G)+ncol(H)),(ncol(G)+1):(ncol(G)+ncol(H))]))
  }
  
  sdAlpha <- matrix(t(G %*% sdPsi),nrow=p0,ncol=r)
  sdBeta  <- matrix(H %*% sdPhi,nrow=p1,ncol=r)
  vSdPi   <- kronecker(beta,diag(p0)) %*% vec(sdAlpha)
  sdPi    <- abs(matrix(vSdPi,nrow=p0,ncol=p1))
  
  covarList           <- list(hessian=hessian)
  covarList$iHessian  <- solve(iHessian)
  covarList$sePi      <- sdPi / sqrt(T)
  covarList$sdAlpha   <- sdAlpha
  covarList$sdBeta    <- sdBeta
  covarList$seAlpha   <- sdAlpha / sqrt(T)
  covarList$seBeta    <- sdBeta / sqrt(T)
  
  return(covarList)
}
