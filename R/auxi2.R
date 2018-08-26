auxI2 <- function(ft, r) {
  p0 <- ncol(ft$R0)
  p1 <- ncol(ft$R1)
  Time <- nrow(ft$R0)
  s1 <- if (is.null(ft$tau)) 0
  else ncol(ft$tau) - r
  s2 <- p0 - r - s1
  ft$beta <- ft$tau %*% ft$rho
  kSi <- -crossprod(ft$kappa, Null(ft$rho))
  eTa	<- crossprod(Null(ft$beta), ft$tau %*% Null(ft$rho))
  ft$delta <-	tcrossprod(Null(ft$tau)) %*% ft$psi
  tmpfit <- lm.fit(cbind(ft$R2 %*% ft$beta + ft$R1 %*% ft$delta, ft$R1 %*% ft$tau), ft$R0)
  tmpcoef <- t(matrix(tmpfit$coef, 2 * r + s1, p0))
  ft$alpha <-	tmpcoef[, seq_len(r), drop = FALSE]
  ft$zeta <- tmpcoef[, seq(r + 1, 2 * r + s1, len = r + s1), drop = FALSE]
  
  if (!is.null(ncol(ft$Z3))) {
    ft$Psi <- t(lm.fit(ft$Z3, ft$Z0 - tcrossprod(ft$Z2 %*% ft$beta +
                                                   ft$Z1 %*% ft$delta, ft$alpha) -
                         ft$Z1 %*% tcrossprod(ft$tau, ft$zeta))$coef)
    colnames(ft$Psi) <- colnames(ft$Z3)
    ft$fitted.values <- with(ft, xts(tcrossprod(Z2 %*% beta + Z1 %*% delta, alpha) +
                                       Z1 %*% tcrossprod(tau, zeta) +
                                       tcrossprod(Z3, Psi),
                                     index(Z0)
    )
    )
  }
  else ft$fitted.values <- with(ft, xts(tcrossprod(Z2 %*% beta + Z1 %*% delta, alpha) +
                                          ft$Z1 %*% tcrossprod(ft$tau, ft$zeta),
                                        index(Z0)
  )
  )
  
  if (s1 != 0) {
    ft$beta1 <- Null(ft$beta) %*% eTa
    ft$alpha1 <- Null(ft$alpha) %*% kSi
    rownames(ft$beta1) <- colnames(ft$R2)
    colnames(ft$beta1) <- paste('beta1(', 1:ncol(ft$beta1), ')', sep = '')
    rownames(ft$alpha1)	<- colnames(ft$R0)
    colnames(ft$alpha1) <- paste('alpha1(', 1:ncol(ft$alpha1), ')', sep = '')
  }
  if (s2 != 0) {
    ft$beta2 <- Null(ft$beta) %*% Null(eTa)
    ft$alpha2 <- Null(ft$alpha) %*% Null(kSi)
    rownames(ft$beta2) <- colnames(ft$R2)
    colnames(ft$beta2) <- paste('beta2(', 1:ncol(ft$beta2), ')', sep = '')
    rownames(ft$alpha2) <- colnames(ft$R0)
    colnames(ft$alpha2) <- paste('alpha2(', 1:ncol(ft$alpha2), ')', sep = '')
  }
  
  if (r > 0) {
    colnames(ft$beta) <- paste('beta(', 1:ncol(ft$beta), ')', sep = '')
    colnames(ft$alpha) <- paste('alpha(', 1:ncol(ft$alpha), ')', sep = '')
    rownames(ft$beta) <- c(colnames(ft$R2))
  }
  
  ft$residuals <- ft$Z0 - ft$fitted.values
  colnames(ft$residuals) <- colnames(ft$Z0)
  class(ft$residuals) <- c('I2res', class(ft$residuals))
  
  class(ft) <- c('I2', 'I1')
  
  return(ft)
}
