anova.I1 <- function(object, ..., df) {
  if(length(objects <- list(object, ...)) == 1) stop('Nothing to test against')
  test <- 2 * (logLik(objects[[1]]) - logLik(objects[[2]]))
  return(data.frame(Value = test,
                    Test = paste('ChiSq(', df, ')', sep = ''),
                    p.value = 1 - pchisq(test, df)
  )
  )
}

