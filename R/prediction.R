prediction.is <- function(obj, data, exogenous, ...) {
  # Function that creates the predicted values of the model 
  # given in obj using the data provided in data
  #
  # Args
  #   obj: an object of type I1
  #   xtsData: an xts object with the data to use for prediction 
  #
  # Returns
  #   predition.I1 type object 
  #
  
  # Central figures
  p   <- ncol(data)
  k   <- obj$call$lags
  n   <- nrow(obj$Z2)
  EOS <- tail(index(obj$Z0),1)  

  #-------------------------------------
  # Construct out of sample prediction
  #-------------------------------------
  # Replicate call to CIModelFrame to construct Z0,Z1,Z2
  mf           <- obj$call
  mf[['data']] <- data
  if(!is.null(exogenous)) { mf[['exogenous']] <- exogenous }
  mf[[1L]]     <- as.name("CIModelFrame")
  mf           <- eval(mf, parent.frame())
  
  pZ0 <- Z0  <- mf$Z0
  pZ1 <- Z1  <- mf$Z1
  pZ2 <- Z2  <- mf$Z2

  for(i in 1:nrow(Z0)) {
    if( i <= n ) {
      if(!is.null(obj$Z2)) {
        pZ1[i,] <- Z1[i,]
        pZ2[i,] <- Z2[i,]
        pZ0[i,] <- pZ1[i,] %*% tcrossprod(obj$beta,obj$alpha) + tcrossprod(pZ2[i,],obj$Psi)
      }
      else {
        pZ1[i,] <- Z1[i,]
        pZ0[i,] <- pZ1[i,] %*% tcrossprod(obj$beta,obj$alpha)
      }
    }
    else {
      if(!is.null(obj$Z2)) {
        pZ1[i,] <- cbind(coredata(pZ1[i-1,1:p] + pZ0[(i-1),]),coredata(Z1[i,(p+1):ncol(Z1)]))
        pZ2[i,] <- cbind(coredata(pZ0[i-1,]),coredata(pZ2[i-1,1:p*(k-2)]),coredata(Z2[i,(p*(k-1)+1):ncol(Z2)]))
        pZ0[i,] <- pZ1[i,] %*% tcrossprod(obj$beta,obj$alpha) + tcrossprod(pZ2[i,],obj$Psi)
      } else {
        pZ1[i,] <- cbind(coredata(pZ1[i-1,1:p] + pZ0[i-1,]),coredata(Z1[i,(p+1):ncol(Z1)]))
        pZ0[i,] <- pZ1[i,] %*% tcrossprod(obj$beta,obj$alpha)
      }
    }
  }
  
  pX <- pZ0+pZ1[,1:p]
  
  # The full data series
  Z0 <- mf$Z0
  X  <- mf$X
  
  predObj <- list(PredX=pX,
                  ActuX=X,
                  EndOfSample=EOS)
  
  class(predObj) <- "prediction.I1"
  
  return(predObj) 
}

prediction.oos <- function(obj, data, steps) {
  # function that generates out of sample predictions
  # a specified number of steps ahead out of the given sample
  # in data.
  #
  # Args:
  #   data: an xts obj containing the data
  #   steps: the number of steps ahead.
  #
  # Returns:
  #   predList: a prediction.I1 object with predicted    
  #             and actual values
  #
  
  # Central figures
  p   <- ncol(data)
  k   <- obj$call$lags
  n   <- nrow(obj$Z2)
  BOS <- head(index(obj$Z0),1)
  EOS <- tail(index(obj$Z0),1)  
  
  # Extend the data with n periods of NA's
  timeIndex <- index(data)
  per       <- periodicity(data)

  
  
  #-------------------------------------
  # Construct out of sample prediction
  #-------------------------------------
  # Replicate call to CIModelFrame to construct Z0,Z1,Z2
  mf           <- obj$call
  mf[['data']] <- data
  mf[[1L]]     <- as.name("CIModelFrame")
  mf           <- eval(mf, parent.frame())
  
  pZ0 <- Z0  <- mf$Z0
  pZ1 <- Z1  <- mf$Z1
  pZ2 <- Z2  <- mf$Z2
  
  for(i in 1:nrow(Z0)) {
    if( i <= n ) {
      if(!is.null(obj$Z2)) {
        pZ1[i,] <- Z1[i,]
        pZ2[i,] <- Z2[i,]
        pZ0[i,] <- pZ1[i,] %*% tcrossprod(obj$beta,obj$alpha) + tcrossprod(pZ2[i,],obj$Psi)
      }
      else {
        pZ1[i,] <- Z1[i,]
        pZ0[i,] <- pZ1[i,] %*% tcrossprod(obj$beta,obj$alpha)
      }
    }
    else {
      if(!is.null(obj$Z2)) {
        pZ1[i,] <- cbind(coredata(pZ1[i-1] + pZ0[(i-1),]),coredata(Z1[i,(p+1):ncol(Z1)]))
        pZ2[i,] <- cbind(coredata(pZ0[i-1,]),coredata(pZ2[i-1,1:p*(k-2)]),coredata(Z2[i,(p*(k-1)+1):ncol(Z2)]))
        pZ0[i,] <- pZ1[i,] %*% tcrossprod(obj$beta,obj$alpha) + tcrossprod(pZ2[i,],obj$Psi)
      } else {
        pZ1[i,] <- cbind(coredata(pZ1[i-1,1:p] + pZ0[i-1,]),coredata(Z1[i,(p+1):ncol(Z1)]))
        pZ0[i,] <- pZ1[i,] %*% tcrossprod(obj$beta,obj$alpha)
      }
    }
  }
  
  pX <- pZ0+pZ1[,1:p]
  
  # The full data series
  Z0 <- mf$Z0
  X  <- mf$X
  
  predObj <- list(PredX=pX,
                  ActuX=X,
                  EndOfSample=EOS)
  
  class(predObj) <- "prediction.I1"
  
  return(predObj) 
}