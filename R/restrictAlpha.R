restrictAlpha <- function(obj, A.matrix, type = 'restr') {
  if(!inherits(A.matrix, 'matrix')) stop('Restriction matrix must be a matrix')
  else if (type == 'restr') {
    H.bar <- A.matrix %*% solve(crossprod(A.matrix))
    H.ort <- Null(A.matrix)
    tmp.rrr <- rrr(obj$R0 %*% H.bar, obj$R1, obj$R0 %*% H.ort, ncol(obj$beta))
    obj$H.alpha <- A.matrix
    obj$alpha <- A.matrix %*% tmp.rrr$alpha
    obj$beta <- tmp.rrr$beta
  }
  else if (type == 'known') {
    a.bar <- A.matrix %*% solve(crossprod(A.matrix))
    a.ort <- Null(A.matrix)
    a.bar.ort <- Null(a.bar)
    if (ncol(A.matrix) == ncol(obj$alpha)) {
      obj$alpha <- A.matrix
      obj$beta <- solve(crossprod(obj$R1), crossprod(obj$R1, obj$R0 %*% a.bar))
    }
    else {
      tmp.rrr <- rrr(obj$R0 %*% a.bar.ort,
                     obj$R1,
                     ,
                     ncol(obj$beta) - ncol(A.matrix)
      )
      tmp.X <- cbind(obj$R1,
                     obj$R0 %*%
                       a.bar.ort - obj$R1 %*%
                       tcrossprod(tmp.rrr$beta, tmp.rrr$alpha)
      )
      tmp.fit <- solve(crossprod(tmp.X), crossprod(tmp.X, obj$R0 %*% a.bar))
      obj$alpha <- cbind(A.matrix, a.ort %*% tmp.rrr$alpha)
      obj$beta <- cbind(tmp.fit[1:ncol(obj$Z1), ], tmp.rrr$beta)
    }
  }
  
  obj <- auxI1(obj)
  
  return(obj)
}
