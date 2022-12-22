
delta_coef <- function(formula, data, start, end) {
  lm(formula, data = data[start:end,])$coeff
}

vcovs <- function(formula, data, start, end) {
  sandwich::kernHAC(lm(formula, data[start:end,]))
}


#' Calculate F-statistic value
#' @param formula y~x formula object
#' @param data cci dataframe
#' @param breaks number of breaks to test for
#' @param ssr SSR Dataframe
fstats <- function(formula, data, breaks, ssr) {

  q <- ncol(model.matrix(formula, data = data))

  obs <- nrow(data)

  ssr_result <- min.ssr(breaks, ssr)

  deltas <- unlist(purrr::map2(ssr_result$start, ssr_result$end, ~delta_coef(formula, data, .x, .y)))
  mdelta <- as.matrix(deltas)

  covs <- purrr::map2(ssr_result$start, ssr_result$end, ~vcovs(formula, data, .x, .y))
  mvcov <- diag(covs)

  list(
    fstat = 1/obs * ((obs-(breaks+1)*q)/breaks*q) * t(mdelta) %*% chol2inv(chol(mvcov)) %*% mdelta,
    breaks = breaks,
    ssr = ssr_result$ssr,
    q = q
  )

}
