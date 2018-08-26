testLagLength <- function(fit, lag.max = 8) {
  
  if (class(fit) != 'I1') stop('Input must have class I1')
  
  ##	m <- match(c('data', 'lags'), names(fit$call))
  Time <- nrow(fit$X) - lag.max
  
  loglik <- Regr <- rep(NA, lag.max)
  infocr <- matrix(NA, lag.max, 2)
  
  for (i in lag.max:1) {
    fit$call[['data']] <- if (i == lag.max) fit$X else fit$X[-seq_len(lag.max - i),]
    fit$call[['lags']] <- i
    tmp <- eval(fit$call, envir = attr(fit$call, 'environment'))
    loglik[i] <- logLik(tmp)
    Regr[i] <- ncol(tmp$Z1) + ifelse(is.null(tmp$Z2), 0, ncol(tmp$Z2))
    infocr[i,] <- infocrit(tmp)[1:2]
  }
  
  z <- list()
  
  z$model.sum	<- data.frame(model = paste('VAR(', lag.max:1, ')', sep = ''),
                            k	= lag.max:1,
                            t	= Time,
                            Regr = rev(Regr),
                            loglik = formatC(rev(loglik), digits = 3, format = 'f'),
                            SC = formatC(rev(infocr[,1]), digits = 3, format = 'f'),
                            HQ = formatC(rev(infocr[,2]), digits = 3, format = 'f')
  )
  
  lagtest <- data.frame()
  for (i in lag.max:1) {
    for (j in (lag.max - 1):1) {
      if (i > j) {
        test <- sprintf('VAR(%i) << VAR(%i)', j, i)
        typetmp	<- ncol(tmp$X) * (Regr[i] - Regr[j])
        type <- sprintf('ChiSqr(%i)', typetmp)
        value <- 2 * (loglik[i] - loglik[j])
        pvalue <- 1 - pchisq(value, typetmp)
        
        z$lagtest <- rbind(z$lagtest,
                           data.frame(test,
                                      type,
                                      value = formatC(value, digits = 3, format = 'f'),
                                      p.value = formatC(pvalue, digits = 3, format = 'f')
                           )
        )
      }
    }
  }
  
  ##    return(lapply(z, format, digits = 6, scientific = F))
  return(z)
}
