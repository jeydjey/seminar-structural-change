.datatable.aware = TRUE


#' calculate a parallel sum
#' @param ... vectors to summarise
#' @param na.rm na.rm
#' @export
psum <- function(..., na.rm = FALSE) {
  dots <- rlang::dots_list(...)

  m <- do.call(cbind, dots)

  rowSums(m, na.rm = na.rm)

}

#' Calculate \overline(Z) matrix
#'
#' @param z regressors as matrix
#' @param bpos vector of positions where break occurs after
#' @export
zbar <- function(z, bpos) {

  r <- dim(z)[1]
  q <- dim(z)[2]

  s <- length(bpos)+1
  b <- c(1, bpos+1, r+1)
  cs <- seq(1, s*q+1, by = q)

  out <- matrix(0, r, s*q)

  for(i in 1:s) {

    ri <- b[i]:(b[i+1]-1)
    ci <- cs[i]:(cs[i+1]-1)
    out[ri, ci] <- z[ri,,drop = F]

  }

  return(out)

}
