tracemckinnon <- function(npr, type = c('trace', 'lambda'), dettrend = c('none', 'mean', 'drift'), tlev, tst) {
  ## Function is just an interface to McKinnons Fortran function. Remember to ask
  ## for accept of use.
  
  npr <- if (is.vector(npr) & all(npr %in% 0:12)) as.integer(npr)
  else stop('npr must be a vector of integer values between 0 and 12')
  itt <- if (type == 'lambda') as.integer(1)
  else if (type == 'trace') as.integer(2)
  else stop('Test type must be either trace or lambda')
  itv <- if (dettrend == 'none') as.integer(0)
  else if (dettrend == 'mean') as.integer(3)
  else if (dettrend == 'drift') as.integer(4)
  else stop('Deterministic specification must be either none, mean or trend')
  if (!missing(tlev)) {
    arg <- if (all(0 < tlev) & all(tlev < 0.5)) as.numeric(tlev)
    else stop('Testlevel must be between 0.0001 and 0.5')
    nc <- as.integer(1)
  }
  else if (!missing(tst)) {
    arg <- if (all(0 < tst) & length(tst) == length(npr)) as.numeric(tst)
    else stop('All teststatistic must be positive and in same number as restrictions')
    nc <- as.integer(2)
  }
  val <- numeric(1)
  rpath <- system.file(package='civecm')
  
  if (nc == 1) {
    ans <- matrix(NA, length(npr), length(arg))
    colnames(ans) <- paste('Cr', substr(arg, 2, 5), sep = '')
    for (i in seq_along(npr)) {
      for (j in seq_along(arg)) {
        ans[i, j] <- .Fortran('johval',
                              npr[i],
                              itt,
                              itv,
                              nc,
                              isave = integer(1),
                              arg[j],
                              val = val,
                              rpath = rpath,
                              as.integer(nchar(rpath))
        )$val
      }
    }
  }
  else if (nc == 2) {
    ans <- matrix(NA, length(npr), 1)
    colnames(ans) <- 'a.p.value'
    for (i in seq_along(npr)) {
      ans[i, 1] <- .Fortran('johval',
                            npr[i],
                            itt,
                            itv,
                            nc,
                            isave = integer(1),
                            arg[i],
                            val = val,
                            rpath = rpath,
                            as.integer(nchar(rpath))
      )$val
    }
  }
  return(ans)
}
