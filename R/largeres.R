largeres <- function(fit, abs.value) {
  res <- rstandard(fit)
  
  if (missing(abs.value)) idx <- cbind(apply(abs(res), 2, which.max), seq(ncol(res)))
  else idx <- which(abs(res) > abs.value, T)
  dates <- index(res)[idx[,1]]
  variable <- apply(idx, 1, function(x) colnames(res)[x[2]])
  values <- as.matrix(res)[idx]
  
  return(data.frame(dates = dates, var = variable, values = values))
}
