boot <- function(bootfct, obj, reps, wild = TRUE, cluster = NULL) {
  reslen <- nrow(obj$residuals)
  sf <- function(dummy) {
    newcall	<- obj$call
    if (wild)
      newcall[['data']] <- simulate(obj, res = obj$residuals * rnorm(reslen))
    else
      newcall[['data']] <- simulate(obj,
                                    res = obj$residuals[sample(seq(reslen),
                                                               replace = TRUE),]
      )
    return(bootfct(eval(newcall, envir = attr(obj$call, 'environment'))))
  }
  if (is.null(cluster)) {
    if (suppressWarnings(require('pbapply'))) ans <- pblapply(seq(reps), sf)
    else ans <- lapply(seq(reps), sf)
  }
  else {
    clusterEvalQ(cluster, library(civecm))
    if (wild) clusterSetupRNG(cluster, type = 'RNGstream')
    ans <- parLapply(cluster, seq(reps), sf)
  }
  return(ans)
}