#' @title A function to print matrices of coefficient estimates with se og t.values

coefmatrix <- function(coef, stat1, digits = 3, sig.dig = c(TRUE, FALSE)) {
  m.coef <- as.matrix(coef)
  m.stat1 <- as.matrix(stat1)
  mchar <- max(nchar(cbind(m.coef, m.stat1)))
  if (any(dim(m.coef)!=dim(m.stat1))) stop('Coef and stat matrices must have same dimension')
  if (length(sig.dig) != 2)
    if (length(sig.dig) == 1) sig.dig <- rep(sig.dig, 2)
  else stop('Argument sig.dig must have length two')
  m.out <- matrix(, 0, ncol(coef))
  for (i in seq_len(nrow(coef))) {
    m.out <- rbind(m.out,
                   formatC(m.coef[i, ],
                           width = -1,
                           digits = digits,
                           format = ifelse(sig.dig[1], 'fg', 'f')
                   ),
                   paste('(',
                         formatC(m.stat1[i, ],
                                 width = -1,
                                 digits = digits,
                                 format = ifelse(sig.dig[2], 'fg', 'f')),
                         ')',
                         sep = '')
    )
  }
  maxi <- apply(m.out, 2, nchar)
  li <- apply(m.out, 2, function(x) regexpr('.', x, fixed = TRUE))
  ri <- maxi - li
  maxcolr <- apply(ri, 2, max)
  maxcoll <- apply(li, 2, max)
  for (i in seq_len(nrow(m.out))) {
    for (j in seq_len(ncol(m.out))) {
      m.out[i, j]	<- paste(paste(rep(' ', maxcoll[j] - li[i, j]),
                                 collapse = ''
      ),
      m.out[i, j],
      paste(rep(' ', maxcolr[j] - ri[i, j]),
            collapse = ''
      ),
      sep = ''
      )
    }
  }
  if (!is.null(rownames(coef))) rownames(m.out) <- c(rbind(rownames(coef), ""))
  colnames(m.out) <- colnames(coef)
  return(m.out)
}
