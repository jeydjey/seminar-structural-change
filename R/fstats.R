
delta_coef <- function(formula, data, from, to) {
  lm(formula, data = data[from:to,])$coeff
}

vcovs <- function(formula, data, from, to) {
  sandwich::kernHAC(lm(formula, data[from:to,]))
}

fcrit <- function(sig, q, breaks) {

  #TODO load crit values and return as requested

}

dcrit <- function(sig, q, type = c("ud", "wd")) {

  #TODO load crit values and return as requested

}


#' Calculate F-statistic value
#' @param formula y~x formula object
#' @param data cci dataframe
#' @param breaks number of breaks to test for
#' @param ssr SSR Dataframe
#' @export
fstats <- function(formula, data, breaks, ssr) {

  z <- model.matrix(formula, data = data)
  q <- ncol(z)

  obs <- nrow(data)

  ssr_result <- min_ssr(breaks, ssr)

  deltas <- unlist(purrr::map2(ssr_result$from, ssr_result$to, ~delta_coef(formula, data, .x, .y)))
  mdelta <- as.matrix(deltas)

  mr <- matrix(0, breaks, breaks+1)
  for(i in 1:breaks) {
    mr[i,i] <- -1
    mr[i,i+1] <- 1
  }
  mr <- mr %x% diag(q)

  covs <- purrr::map2(ssr_result$from, ssr_result$to, ~vcovs(formula, data, .x, .y))
  mvcov <- diag(covs)

  f <- 1/obs * ((obs-(breaks+1)*q)/breaks*q) * t(mdelta) %*% t(mr) %*% solve(mr %*% mvcov %*% t(mr)) %*% mr %*% mdelta

  list(
    fstat = f[1,1],
    breaks = breaks,
    ssr = ssr_result$ssr,
    q = q
  )

}

#' Double maximum test for up to the number of breaks
#' @param formula y~x formula object
#' @param data cci dataframe
#' @param breaks number of 1:breaks to test for
#' @param ssr SSR Dataframe
#' @param sig significance level
#' @export
dmax_stats <- function(formula, data, breaks, ssr, sig = "5%") {

  fstat <- unlist(purrr::map(1:breaks, ~fstats(formula, data, .x, ssr)$fstat))

  ud <- max(fstat)

  #TODO wd implementing and crit values

}
