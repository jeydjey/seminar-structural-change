
delta_coef <- function(formula, data, from, to) {
  lm(formula, data = data[from:to,])$coeff
}

vcovs <- function(formula, data, from, to) {
  sandwich::kernHAC(lm(formula, data[from:to,]), kernel = "Quadratic Spectral", approx = "AR(1)")
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
#' @export
fstats <- function(formula, data, breaks, ssr, conf = 0.95, trim = 0.1) {

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
  mdelta <- as.matrix(deltas)

  mr <- matrix(0, breaks, breaks+1)
  for(i in 1:breaks) {
    mr[i,i] <- -1
    mr[i,i+1] <- 1
  }
  mr <- mr %x% diag(q)

  covs <- purrr::map2(ssr_result$from, ssr_result$to, ~vcovs(formula, data, .x, .y))
  mvcov <- diag(covs)

  f <- 1/obs * (((obs-(breaks+1)*q)/breaks*q) %*% t(mdelta) %*% t(mr) %*% solve(mr %*% mvcov %*% t(mr)) %*% mr %*% mdelta)

  crit <- fcrit(conf = 0.95, trim = trim, q = q, k = breaks)

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
#' @export
#' @importFrom magrittr %>% %$%
dmax_stats <- function(formula, data, breaks, ssr, conf = 0.95, trim = 0.1) {

  z <- model.matrix(formula, data = data)
  q <- ncol(z)

  fstat <- dplyr::bind_rows(purrr::map(1:breaks, ~fstats(formula, data, .x, ssr) %$%
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
#' @export
#' @importFrom magrittr %>%
lstats <- function(formula, data, l, ssr, conf = 0.95, trim = 0.1) {

  y <- model.response(model.frame(formula, data = data))
  q <- ncol(z)

  ssr_result <- min_ssr(l, ssr)

  fstat <- dplyr::bind_rows(
    purrr::map2(ssr_result$from, ssr_result$to, function(.x, .y) {

      ssr_2 <- ssr %>%
        dplyr::filter(to <= .y, to_from >= (.x-1)) %>%
        dplyr::mutate(to = to - (.x-1),
                      to_from = to_from - (.x-1))

      fstats(formula, data[.x:.y,], 1, ssr_2, conf, trim)


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
#' @export
#' @importFrom magrittr %>%
seq_test <- function(formula, data, m, ssr, conf = 0.95, trim = 0.1) {

  y <- model.response(model.frame(formula, data = data))
  q <- ncol(z)

  result <- min_ssr(0, ssr)

  for(i in 1:m) {

    ltest <- lstats(formula, data, i-1, ssr, conf, trim)

    if(!ltest$significant)
      return(result)

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

        stats <- fstats(formula, data[.y$from:.y$to,], 1, ssr_2, conf, trim)

        list(
          id = .y$id,
          fstat = stats$fstat,
          breaks_from = stats$breaks_from,
          breaks_to = stats$breaks_to,
          ssr = stats$ssr
        )

      }) %>% dplyr::bind_rows() %>%
      dplyr::slice_max(fstat)

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

