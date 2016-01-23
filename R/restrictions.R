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

restrictBeta <- function(obj, H.matrix) {
  # Reestimates the model under the restrictions provided in H.matrix 
  # This function reestimates alpha implicitly under the restrictions on
  # on beta and hence will annihilate any restrictions previously imposed 
  # on alpha. If restrictions on both alpha and beta are needed, please
  # use the function restrictLongRun.
  #
  # Args: 
  #   obj: an object of type I1
  #   H.matrix: a matrix with restrictions to impose on beta.
  #     (note that this restriction matrix is not equivalent to those
  #      used in restrictLongRun which impose restrictions on vec(beta)   
  #      rather that beta itself)
  #
  # Returns:
  #   obj: an object of type I1 with the restrictions imposed on beta.
  #
  
  if (length(H.matrix) == 1) H.matrix	<- H.matrix[[1]]
  if (inherits(H.matrix, 'matrix')) {
    temp.rrr   <- rrr(obj$Z0, obj$Z1 %*% H.matrix, obj$Z2, ncol(obj$beta))
    obj$H.beta <- H.matrix
    obj$beta   <- H.matrix %*% temp.rrr$beta
    ## obj$beta <- betaNorm(obj$beta, obj$S11)
    obj$alpha  <- t(lm.fit(obj$R1%*%obj$beta, obj$R0)$coef)
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
    
    i       <- 1
    loglik0 <- logLik(obj)
    iter    <- 1
    while(1) {
      temp.rrr <- rrr(obj$R0,
                      obj$R1 %*% H.matrix[[i]],
                      obj$R1 %*% obj$beta[, -i, drop = F],
                      1
      )
      
      obj$beta[, i] <- H.matrix[[i]] %*% temp.rrr$beta
      
      if (i == ncol(obj$beta)) {
        ##obj$beta		<- betaNorm(obj$beta, obj$S11)
        tmp.lm        <- lm.fit(obj$R1%*%obj$beta, obj$R0)
        obj$alpha     <- t(tmp.lm$coef)
        obj$residuals <- tmp.lm$residuals
        
        loglik1 <- logLik(obj)
        if (loglik1<loglik0 && iter != 1) {
          print(c(loglik0, loglik1))
          stop('Loglik is decreasing. Try to impose identifying restrictions')
        }
        if (abs((loglik0 - loglik1) / loglik0) < 1e-6) {
          cat('Convergence in', iter, 'iterations.\n')
          break
        }
        
        loglik0	<- loglik1
        i       <- 1
        iter    <- iter + 1
      }
      else i <- i + 1
    }
  }
  obj	<- auxI1(obj)
  
  return(obj)
}

restrictLongRun <- function(obj,Hb,hb,Ha,beta_i,alpha_i) {
  # Function that restricts both subelements of the long run matrix
  # The function applies the general Boswijk/Doornik principle for
  # imposing linear restrictions on the long run parameters.
  # See Boswijk and Doornik (2004) : "Identifying, Estimating and
  # Testing restricted cointegration systems"
  #
  # Args
  #   obj: an I1 object with restricted parameter estimates
  #   Hb: matrix with linear restrictions on vec(beta)
  #   hb: a vector of corresponding fixed parameter values
  #   Ha: matrix with linear restrictions on vec(alpha)
  #   alpha_i: a matrix with initial values for alpha
  #   beta_i: a matrix with initial values for beta
  #
  # Returns
  #   obj: an object of type I1  
  #
  
  p0               <- ncol(obj$Z0)
  p1               <- ncol(obj$Z1)
  p2               <- ifelse(is.null(obj$Z2), 0, ncol(obj$Z2))
  p                <- p0
  r                <- ncol(alpha_i)
  p1               <- nrow(beta_i)
  obj_i            <- obj
  obj_i$alpha      <- alpha_i
  obj_i$beta       <- beta_i
  obj_i            <- auxI1(obj_i)
  obj_i$Hb         <- Hb
  obj_i$hb         <- hb
  obj_i$Ha         <- Ha
  obj_i$LongRunRes <- TRUE
  S00              <- crossprod(obj_i$R0,obj_i$R0) / nrow(obj_i$R0)   
  S01              <- crossprod(obj_i$R0,obj_i$R1) / nrow(obj_i$R0)
  S11              <- crossprod(obj_i$R1,obj_i$R1) / nrow(obj_i$R1)
  alpha_prime_i    <- t(alpha_i)  
  PI_hat_LS        <- solve(S11,t(S01))
  Omega_i          <- (S00 - S01 %*% tcrossprod(beta_i,alpha_i) 
                       - t(S01 %*% tcrossprod(beta_i,alpha_i))
                       + tcrossprod(alpha_i,beta_i) %*% S11 %*% tcrossprod(beta_i,alpha_i))
  invOmega_i <- solve(Omega_i)
  
  loglik_i <- logLik(obj_i)
  loglik_1 <- loglik_i - 1
  
  i <- 1
  while( (loglik_i - loglik_1) > 1e-6 && i < 10000) {
    # Reset loglikelihood
    loglik_1 <- loglik_i
    
    # Estimate alpha
    psi_i          <- solve(crossprod(Ha, kronecker(invOmega_i,t(beta_i) %*% S11 %*% beta_i) %*% Ha), 
                            crossprod(Ha, kronecker(invOmega_i,t(beta_i) %*% S11))) %*% vec(PI_hat_LS)
    valpha_prime_i <- Ha %*% psi_i
    alpha_prime_i  <- matrix(valpha_prime_i,ncol=p)
    alpha_i        <- t(alpha_prime_i)
    
    # Estimate beta
    if(!is.null(Hb)) {
      phi_i   <- solve(crossprod(Hb, kronecker(crossprod(alpha_i,invOmega_i %*% alpha_i),S11) %*% Hb),
                       crossprod(Hb, kronecker(crossprod(alpha_i,invOmega_i),S11))) %*% (vec(PI_hat_LS) - kronecker(alpha_i,diag(p1)) %*% hb)
      vbeta_i <- Hb %*% phi_i + hb
      beta_i  <- matrix(vbeta_i,ncol=r)
    }
    # Estimate Omega
    Omega_i    <- (S00 - S01 %*% tcrossprod(beta_i,alpha_i) - t(S01 %*% tcrossprod(beta_i,alpha_i))
                   + tcrossprod(alpha_i,beta_i) %*% S11 %*% tcrossprod(beta_i,alpha_i))
    invOmega_i <- solve(Omega_i)
    
    # Calculate likelihood
    obj_i$alpha <- alpha_i
    obj_i$beta  <- beta_i
    
    obj_i    <- auxI1(obj_i)
    loglik_i <- logLik(obj_i)
    i        <- i + 1
  }
  obj <- auxI1(obj_i);
  return(obj)
}

restrictTau <- function(obj, Hmat) {
  if (!inherits(obj, 'I2')) stop('Object must have class I2')
  if (nrow(Hmat) != ncol(obj$R1))
    stop('Restriction matrix has wrong number of rows')
  if (!(ncol(Hmat) %in% seq(ncol(obj$R1))))
    stop('Restriction matrix has wrong number of columns')
  Time <- nrow(obj$R0)
  p0 <- ncol(obj$R0)
  p1 <- ncol(obj$R1)
  s1 <- ncol(obj$beta1)
  s2 <- ncol(obj$beta)
  
  if (p0 - s1 - s2 > 0)
    ans <- .Fortran('tauswitch',
                    as.double(obj$R0),
                    as.double(obj$R1),
                    as.double(obj$R2),
                    as.double(Hmat),
                    time = as.integer(Time),
                    as.integer(p0),
                    as.integer(p1),
                    as.integer(p0 - s1 - s2),
                    as.integer(s2),
                    as.integer(ncol(Hmat)),
                    tau = obj$tau,
                    rho = matrix(0.0, s2 + s1, p0 - s1 - s2),
                    psi = matrix(0.0, p1, s2),
                    kappa = matrix(0.0, s2 + s1, p0 - s2),
                    as.integer(10000), 1e-9, integer(1))
  else
    ans <- restrictBeta(update(obj, s1), Hmat)
  
  obj$tau   <- ans$tau
  obj$rho   <- ans$rho
  obj$psi   <- ans$psi
  obj$kappa <- ans$kappa
  obj	      <- aux.I2(obj, p0 - s1 - s2)
  return(obj)
}
