restrictBeta <- function(obj, H.matrix) {
  if (length(H.matrix) == 1) H.matrix	<- H.matrix[[1]]
  if (inherits(H.matrix, 'matrix')) {
    temp.rrr <- rrr(obj$Z0, obj$Z1 %*% H.matrix, obj$Z2, ncol(obj$beta))
    obj$H.beta <- H.matrix
    obj$beta <- H.matrix %*% temp.rrr$beta
    ## obj$beta <- betaNorm(obj$beta, obj$S11)
    obj$alpha <- t(lm.fit(obj$R1%*%obj$beta, obj$R0)$coef)
  }
  else if (!inherits(H.matrix, 'list'))
    stop('H.matrix must be either a matrix or a list of matrices.')
  else if (length(H.matrix) != ncol(obj$beta))
    stop('The number of elements in H.matrix must match the cointegration rank')
  else {
    for (i in seq_len(ncol(obj$beta))) {
      obj$beta[, i] <- obj$beta %*%
        geigen(crossprod(obj$beta, H.matrix[[i]]) %*%
                 solve(crossprod(H.matrix[[i]]),
                       crossprod(H.matrix[[i]], obj$beta)
                 ),
               crossprod(obj$beta)
        )$vectors[, 1]
    }
    
    i <- 1
    loglik0 <- logLik(obj)
    iter <- 1
    
    while(1) {
      temp.rrr <- rrr(obj$R0,
                      obj$R1 %*% H.matrix[[i]],
                      obj$R1 %*% obj$beta[, -i, drop = F],
                      1
      )
      
      obj$beta[, i] <- H.matrix[[i]] %*% temp.rrr$beta
      
      if (i == ncol(obj$beta)) {
        ##obj$beta		<- betaNorm(obj$beta, obj$S11)
        tmp.lm <- lm.fit(obj$R1%*%obj$beta, obj$R0)
        obj$alpha <- t(tmp.lm$coef)
        obj$residuals <- tmp.lm$residuals
        loglik1 <- logLik(obj)
        if (loglik1<loglik0 && iter != 1) {
          print(c(loglik0, loglik1))
          stop('Loglik is decreasing. Try to impose identifying restrictions')
        }
        if (abs((loglik0 - loglik1) / loglik0) < 1e-20) {
          cat('Convergence in', iter, 'iterations.\n')
          break
        }
        loglik0	<- loglik1
        i <- 1
        iter <- iter + 1
      }
      else i <- i + 1
    }
  }
  obj	<- auxI1(obj)
  
  return(obj)
}
