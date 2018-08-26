#' @title Function borrowed from MASS but with here is the convention that the 
#'        Null of a zero column matrix is the identity

Null <- function(M) {
  tmp <- qr(M)
  set <- if (tmp$rank == 0) seq_len(ncol(M))
  else -seq_len(tmp$rank)
  ans <- if (ncol(M) == 0) diag(nrow(M))
  else qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
  return(ans)
}
