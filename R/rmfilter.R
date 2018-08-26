rmfilter <- function(x, coef, init) {
  if (any(is.na(x))) stop('Series includes NAs')
  if (any(is.na(coef))) stop('Coefficients includes NAs')
  if (is.null(dim(x))) stop('Input series has no dimensions')
  nl <- dim(x)
  if (!is.array(coef)) stop('coef must be an array')
  cdim <- dim(coef)
  if (nl[2] != cdim[2] | cdim[1] != cdim[2]) stop('coef has wrong dimensions')
  if (missing(init)) init <- matrix(0, cdim[3], nl[2])
  else if (any(is.na(init))) stop('Initial values includes NAs')
  else if (any(dim(init) != c(cdim[3], nl[2]))) stop('init matrix has wrong dimensions')
  storage.mode(x) <- 'double'
  ans	<- .Fortran('rmfilter',
                  x,
                  as.integer(nl),
                  as.double(coef),
                  as.integer(cdim),
                  as.double(init)
  )[[1]]
  return(ans)
}
