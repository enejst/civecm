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
  
  obj$tau <- ans$tau
  obj$rho <- ans$rho
  obj$psi <- ans$psi
  obj$kappa <- ans$kappa
  obj	<- auxI2(obj, p0 - s1 - s2)
  return(obj)
}
