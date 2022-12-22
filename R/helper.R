datatable.aware = TRUE


#' calculate a parallel sum
#' @param ... vectors to summarise
#' @param na.rm
#' @export
psum <- function(..., na.rm = FALSE) {
  dots <- rlang::dots_list(...)

  m <- do.call(cbind, dots)

  rowSums(m, na.rm = na.rm)

}
