update.I1 <- function(object, rank, ...) {
  stopifnot(length(rank) %in% c(1, 2))
  if (length(rank) == 2 && rank[2] == 0) {
    object$call[['rank']] <- rank[1]
    ans	<- eval(object$call, envir = attr(object$call, 'environment'))
  }
  else {
    object[['call']][['rank']] <- rank
    ans <- eval(object[['call']], envir = attr(object$call, 'environment'))
  }
  return(ans)
}
