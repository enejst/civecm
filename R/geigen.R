geigen	<- function(A, B, symmetric, only.values = FALSE)
{
	if(missing(symmetric))
	{
		symmetric	<- (all.equal(A, t(A)) & all.equal(B, t(B)))
	}
	else if(!is.logical(symmetric)) stop('Argument symmetric must be logical')
	else if(!symmetric) stop('Only the symmetric case is handled at this point')
	if(any(dim(A) != dim(B))) stop('A and B must have same dimensions')
	jobz	<- ifelse(only.values, 'N', 'V')
	n		<- ncol(A)
	ans1	<- .Fortran('dsygvd', as.integer(1), jobz, 'U', as.integer(n), A = A, as.integer(n),
			B, as.integer(n), values = rep(0.0, n), rep(0.0, 4*(1 + 6*n + 2*n**2)), 
			as.integer(4*(1 + 6*n + 2*n**2)), as.integer(rep(0, 4*(3 + 5*n))), 
			as.integer(4*(3 + 5*n)), info = as.integer(0))
	if(ans1$info != 0) stop('Error in Fortran calculation')
	ans2	<- list()
	ans2$values	<- rev(ans1$values)
	if(!only.values) ans2$vectors <- ans1$A[, rev(seq(ncol(A)))]
	return(ans2)
}