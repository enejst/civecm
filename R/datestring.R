datestring <- function(dates, freq) {
  ans <- if(freq %in% c(4, 12)) sapply(dates, function(x) paste(floor(x), ifelse(freq == 4, ' Q', ' M'), round((x - floor(x))*freq + 1), sep = ''))
  else dates
  return(ans)
}
