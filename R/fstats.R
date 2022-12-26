
delta_coef <- function(formula, data, from, to) {
  lm(formula, data = data[from:to,])$coeff
}

hac <- function(zb, u, prewhite = 1) {

  r <- dim(zb)[1]
  s <- dim(zb)[2]

  b = matrix(0, nrow = s, ncol = 1)
  mb = matrix(0, nrow = s, ncol = d)
  hac

}

vcovs <- function(formula, data, from, to, autocorrelation = FALSE) {

  z <- model.matrix(formula, data)
  x <- model.response(model.frame(formula, data = data))

  #checking data distributions across segments
  data_dist_eq <- all(unlist(purrr::map(seq_len(length(from)-1), ~ks.test(x[from[.x]:to[.x]], x[from[.x+1]:to[.x+1]])$p.value))>=0.05)
  #checking error distribution across segments
  error_dist_eq <- all(unlist(purrr::map(seq_len(length(from)-1), ~ks.test(residuals(lm(formula, data[from[.x]:to[.x],])),
                                                                           residuals(lm(formula, data[from[.x+1]:to[.x+1],])))$p.value))>=0.05)
  zb <- zbar(z, from[-1]-1)

  if(!autocorrelation & !data_dist_eq & error_dist_eq){
    #No serial correlation, different distributions for the data, identical distributions for the errors
    mdl <- lm(formula, data)
    s_hat <- sigma(mdl)^2
    return(s_hat*solve((1/nrow(data)*t(zb)%*%zb)))
  }

  if(!autocorrelation & !data_dist_eq & !error_dist_eq){
    #No serial correlation, different distribution for the data, different distributions for the errors
    s_hat <- diag(
      unlist(purrr::map2(from, to, function(.x,.y) {
        d <- data[.x:.y,]
        mdl <- lm(formula, d)
        return(1/nrow(d) * sum(residuals(mdl)^2))
      }))
    )
    return(s_hat %x% diag(1, dim(z)[2]) %*% solve(t(zb) %*% zb))
  }

  if(!autocorrelation & data_dist_eq & error_dist_eq){
    #No serial correlation, identical distribution for the data, identical variances for the errors
    mdl <- lm(formula, data)
    s_hat <- sigma(mdl)^2
    l_hat <- diag(c(diff(from), nrow(data)-sum(diff(from))))/nrow(data)
    return(s_hat*solve(l_hat %*% (t(zb)%*%zb)))
  }

  if(!autocorrelation & data_dist_eq & !error_dist_eq){
    #No serial correlation, identical distribution for the data, different distributions for the errors
    s_hat <- diag(
      unlist(purrr::map2(from, to, function(.x,.y) {
        d <- data[.x:.y,]
        mdl <- lm(formula, d)
        return(1/nrow(d) * sum(residuals(mdl)^2))
      }))
    )
    l_hat <- diag(c(diff(from), nrow(data)-sum(diff(from))))/nrow(data)
    return(s_hat %x% diag(1, dim(z)[2]) %*% solve(l_hat %x% (t(zb)%*%zb)))
  }

  if(autocorrelation & !data_dist_eq & !error_dist_eq) {
    #serial correlation, different distributions for data and errors
    deltat <- diag(c(diff(from), nrow(data)-sum(diff(from))))
    hac <- diag(unlist(purrr::map2(from, to, ~sandwich::kernHAC(lm(formula, data[.x:.y,]), prewhite = (.y-.x)>3))))
    #return(deltat %*% solve(t(zb) %*% zb) %*% hac %*% solve(t(zb) %*% zb))
    return(hac)
  }

  if(autocorrelation & error_dist_eq) {
    #serial correlation, same distribution for the errors
    hac <- sandwich::kernHAC(lm(formula, data))
    l_hat <- diag(c(diff(from), nrow(data)-sum(diff(from))))/nrow(data)
    #return(nrow(data) * solve(t(zb) %*% zb) %*% (l_hat %x% hac) %*% solve(t(zb) %*% zb))
    return(hac)
  }

  stop(glue::glue("vcov case for autocorrelation = {autocorrelation}, equal data distribution = {data_dist_eq} and equal error distribution = {error_dist_eq} is not defined"),
      call. = FALSE)

}

fcrit <- function(conf, trim, q, k) {

  dplyr::pull(dplyr::filter(readRDS(system.file("extdata", "crit_supf.rda", package = "seminar.bp.strucchange", mustWork = T)),
                conf == !!conf, trim == !!trim, q == !!q, k == !!k), crit)

}

dcrit <- function(conf, trim, q, dmax = c("udmax", "wdmax")) {

  dplyr::pull(dplyr::filter(readRDS(system.file("extdata", "crit_dmax.rda", package = "seminar.bp.strucchange", mustWork = T)),
                            conf == !!conf, trim == !!trim, q == !!q, dmax == !!dmax), crit)
}

lcrit <- function(conf, trim, q, l) {

  dplyr::pull(dplyr::filter(readRDS(system.file("extdata", "crit_seq.rda", package = "seminar.bp.strucchange", mustWork = T)),
                            conf == !!conf, trim == !!trim, q == !!q, l == !!l), crit)

}

#' Calculate F-statistic value
#' @param formula y~x formula object
#' @param data cci dataframe
#' @param breaks number of breaks to test for
#' @param ssr SSR Dataframe
#' @param conf confidence interval
#' @param trim trimming 0.05 - 0.25
#' @param autocorrelation does the data contain autocorrelation?
#' @export
fstats <- function(formula, data, breaks, ssr, conf = 0.95, trim = 0.1, autocorrelation = FALSE) {

  z <- model.matrix(formula, data = data)
  q <- ncol(z)

  obs <- nrow(data)

  ssr_result <- min_ssr(breaks, ssr, init = min(ssr$to_from))

  if(ssr_result$breaks == 0) {
    return(
      list(
        fstat = 0,
        breaks = 0,
        ssr = ssr_result$global_ssr,
        crit = NA_real_,
        significant = FALSE
      )
    )
  }

  deltas <- unlist(purrr::map2(ssr_result$from, ssr_result$to, ~delta_coef(formula, data, .x, .y)))
  mdelta <- matrix(deltas)

  mr <- matrix(0, breaks, breaks+1)
  for(i in 1:breaks) {
    mr[i,i] <- -1
    mr[i,i+1] <- 1
  }
  mr <- mr %x% diag(q)

  mvcov <- vcovs(formula, data, ssr_result$from, ssr_result$to, autocorrelation = autocorrelation)

  f <- ((obs-(breaks+1)*q)/(breaks*q*obs)) %*% t(mdelta) %*% t(mr) %*% solve(mr %*% mvcov %*% t(mr)) %*% mr %*% mdelta

  crit <- fcrit(conf = conf, trim = trim, q = q, k = breaks)

  list(
    fstat = f[1,1],
    breaks = breaks,
    breaks_from = ssr_result$from,
    breaks_to = ssr_result$to,
    ssr = ssr_result$ssr,
    global_ssr = ssr_result$global_ssr,
    crit = crit,
    significant = f[1,1] > crit
  )

}

#' Double maximum test for up to the number of breaks
#' @param formula y~x formula object
#' @param data cci dataframe
#' @param breaks number of 1:breaks to test for
#' @param ssr SSR Dataframe
#' @param conf confidence interval
#' @param trim trimming 0.05 - 0.25
#' @param autocorrelation does the data contain autocorrelation?
#' @export
#' @importFrom magrittr %>% %$%
dmax_stats <- function(formula, data, breaks, ssr, conf = 0.95, trim = 0.1, autocorrelation = FALSE) {

  z <- model.matrix(formula, data = data)
  q <- ncol(z)

  fstat <- dplyr::bind_rows(purrr::map(1:breaks, ~fstats(formula, data, .x, ssr, autocorrelation = autocorrelation) %$%
                               tibble::tibble(fstat = fstat, crit = crit)))

  ud <- max(fstat$fstat)

  wd <- fstat %>%
    dplyr::mutate(fstat = fstat * dplyr::first(crit)/crit) %>%
    dplyr::pull(fstat) %>%
    max()

  ud_crit <- dcrit(conf = conf, trim = trim, q = q, dmax = "udmax")
  wd_crit <- dcrit(conf = conf, trim = trim, q = q, dmax = "wdmax")

  list(
    ud = ud,
    ud_crit = ud_crit,
    ud_significant = ud > ud_crit,
    wd = wd,
    wd_crit = wd_crit,
    wd_significant = wd > wd_crit
  )

}

#' F-test for l+1|l breaks
#' @param formula y~x formula object
#' @param data cci dataframe
#' @param breaks number of 1:breaks to test for
#' @param ssr SSR Dataframe
#' @param conf confidence interval
#' @param trim trimming 0.05 - 0.25
#' @param autocorrelation does the data contain autocorrelation?
#' @export
#' @importFrom magrittr %>%
lstats <- function(formula, data, l, ssr, conf = 0.95, trim = 0.1, autocorrelation = FALSE) {

  z <- model.matrix(formula, data = data)
  q <- ncol(z)

  ssr_result <- min_ssr(l, ssr)

  fstat <- dplyr::bind_rows(
    purrr::map2(ssr_result$from, ssr_result$to, function(.x, .y) {

      ssr_2 <- ssr %>%
        dplyr::filter(to <= .y, to_from >= (.x-1)) %>%
        dplyr::mutate(to = to - (.x-1),
                      to_from = to_from - (.x-1))

      fstats(formula, data[.x:.y,], 1, ssr_2, conf, trim, autocorrelation = autocorrelation)


    })
  )

  fstat <- max(fstat$fstat)

  crit <- lcrit(conf = conf, trim = trim, q = q, l = l)

  list(
    l = l,
    fstat = fstat,
    crit = crit,
    significant = fstat > crit
  )

}

#' Sequential testing for another break up to m
#' @param formula y~x formula object
#' @param data cci dataframe
#' @param breaks number of 1:breaks to test for
#' @param ssr SSR Dataframe
#' @param conf confidence interval
#' @param trim trimming 0.05 - 0.25
#' @param skip Should the test of significance for skip+1|skip be skipped, -1 for no skip of 1|0
#' @param autocorrelation does the data contain autocorrelation?
#' @export
#' @importFrom magrittr %>%
seq_test <- function(formula, data, m, ssr, conf = 0.95, trim = 0.1, skip = -1, autocorrelation = FALSE) {

  z <- model.matrix(formula, data = data)
  q <- ncol(z)

  result <- min_ssr(0, ssr)

  for(i in 1:m) {

    fstat <- purrr::map(1:length(result$from), function(.x) {

        .y <- list(
          from = result$from[.x],
          to = result$to[.x],
          id = result$id[.x],
          ssr = result$ssr[.x]
        )

        ssr_2 <- ssr %>%
          dplyr::filter(to <= .y$to, to_from >= (.y$from-1)) %>%
          dplyr::mutate(to = to - (.y$from-1),
                        to_from = to_from - (.y$from-1))

        stats <- fstats(formula, data[.y$from:.y$to,], 1, ssr_2, conf, trim, autocorrelation = autocorrelation)

        list(
          id = .y$id,
          fstat = stats$fstat,
          breaks_from = stats$breaks_from + .y$from - 1,
          breaks_to = stats$breaks_to + .y$from - 1,
          ssr = stats$ssr
        )

      }) %>% dplyr::bind_rows() %>%
      dplyr::slice_max(fstat)

    if(fstat$fstat[1] < lcrit(conf, trim, q, i-1) && !(i-1 <= skip))
      return(result)

      old_id <- fstat$id[1]
      pos <- which(result$id==old_id)

      result$breaks <- result$breaks+1
      result$from <- c(result$from[-pos], fstat$breaks_from)
      o <- order(result$from)
      result$to <- c(result$to[-pos], fstat$breaks_to)[o]
      result$ssr <- c(result$ssr[-pos], fstat$ssr)[o]
      result$from <- result$from[o]
      result$id <- seq_len(result$breaks+1)
      result$global_ssr <- sum(result$ssr)

  }

  return(result)


}


