#' Function to decide on valid segments in series
#' @param from begin of segment
#' @param to to of segment
#' @param length length of total observation period
#' @param n_breaks number of breaks in data
#' @param min_length minimum length of segment
#' @export
valid_segm <- function(from, to, length, n_breaks, min_length) {

  from <= to &
    ((from-1) %/% min_length + (length-to) %/% min_length) >= n_breaks &
    (to - from + 1) >= min_length &
    ((from-1) %/% min_length > 0 | from == 1)

}


#' calculate the squared recurive residual in a segment
#' @param formula y~x formula object
#' @param from from of segment
#' @param to to of segment
#' @param data dataframe of observations
#' @export
rec_res <- function(formula, from, to, data) {

  model <- lm(formula, data = data[from:to,])

  list(
    to = from:to,
    res_sq = c(0, strucchange::recresid(model)^2)
  )
}



#' calculate sum of squared residuals
#' @details Calculates all SSR under consideration of the minimal length requirement
#' @param forumla x~y formulate object
#' @param data dataframe of observations
#' @param trim trimming for minimal length, either percentage or absolute value
#' @export
#' @importFrom magrittr %>%
ssr <- function(formula, data, trim = 0.1) {

  ssr_mtx <- matrix(nrow = nrow(data), ncol = nrow(data))

  min_length <- if(trim<1) floor(trim * nrow(data)) else max(trim, floor(0.05 * nrow(data)))

  df <- tibble::tibble(
    from = c(row(ssr_mtx)),
    to = c(col(ssr_mtx)),
    ssr = mapply(valid_segm, row(ssr_mtx), col(ssr_mtx), nrow(data), 0, min_length)
  )

  # calculate residuals necessary for sum squared residuals, avoiding recursive calculation
  dt <- df %>%
    dplyr::filter(ssr) %>%
    dplyr::group_by(from) %>%
    dplyr::slice_max(to) %>%
    dplyr::ungroup() %>%
    data.table::as.data.table()

  dt <- dt[,
           rec_res(formula, from, to, data),
           by = c("from")
  ]

  # SSR for from to to sequence
  residuals <- dt %>%
    tibble::as_tibble() %>%
    dplyr::group_by(from) %>%
    dplyr::mutate(
      ssr_val = cumsum(res_sq)
    ) %>%
    dplyr::ungroup()

  # join SSR with valid segments
  ssr_df <- df %>%
    dplyr::filter(ssr) %>%
    dplyr::left_join(
      residuals %>%
        dplyr::select(from, to, ssr_val),
      by = c("from", "to")
    )

  # create permutations of possible segment combinations
  ssr_df <- ssr_df %>%
    dplyr::mutate(
      to_from = from-1
    )

  return(ssr_df %>%
    dplyr::select(to, to_from, ssr = ssr_val)
  )

}


#' Dynamic programming algorithm for the global minimum of SSR for a number of breaks
#' @param breaks number of breaks
#' @param ssr ssr table
#' @param init begin of first observation - 1
#' @export
#' @importFrom magrittr %>%
min_ssr <- function(breaks, ssr, init = 0) {

  # dataframe of permutations
  permutations <- tibble::tibble(
    to_0 = init
  )

  # join everything together at the breakpoints
  # and filter after one-break for minimal global ssr before next iteration
  # same as recursive solution but not recursive
  for(i in 1:(breaks+1)) {

    permutations <- permutations %>%
      dplyr::left_join(ssr %>%
                         dplyr::rename(setNames(c("to", "ssr"), c(paste0("to_", i), paste0("ssr_", i)))),
                       by = setNames("to_from", paste0("to_", i-1))) %>%
      dplyr::mutate(global_ssr = psum(!!!rlang::syms(paste("ssr", 1:i, sep = "_"))))

    if(i > 2) {
      permutations <- permutations %>%
        dplyr::group_by(!!rlang::sym(paste("to", i, sep = "_"))) %>%
        dplyr::slice_min(global_ssr, n = 1, with_ties = T) %>%
        dplyr::ungroup()
    }

  }

  # filter only permutation across full observation period
  permutations <- permutations %>%
    dplyr::slice_max(!!rlang::sym(paste0("to_", breaks+1))) %>%
    dplyr::slice_min(global_ssr)

  segments <- dplyr::bind_cols(
    permutations %>%
      dplyr::select(dplyr::all_of(paste("to", seq_len(breaks + 1), sep = "_"))) %>%
      tidyr::pivot_longer(cols = everything(), values_to = "to") %>%
      dplyr::select(to),
    permutations %>%
      dplyr::select(dplyr::all_of(paste("ssr", seq_len(breaks + 1), sep = "_"))) %>%
      tidyr::pivot_longer(cols = everything(), values_to = "ssr") %>%
      dplyr::select(ssr)
  ) %>%
    dplyr::mutate(from = dplyr::lag(to + 1, default = 1),
                  id = dplyr::row_number())

  if(nrow(permutations)==0) {

    return(
      list(
        breaks = 0,
        from = 1,
        to = max(ssr$to),
        id = 1,
        ssr = ssr %>% dplyr::filter(to_from == 0, to == max(to)) %>% dplyr::pull(ssr),
        global_ssr = ssr %>% dplyr::filter(to_from == 0, to == max(to)) %>% dplyr::pull(ssr)
      )
    )
  }

  list(
    breaks = breaks,
    from = segments$from,
    to = segments$to,
    id = segments$id,
    ssr = segments$ssr,
    global_ssr = sum(segments$ssr)
  )
}

#' Calculate the Bayesian Information Criteron for 0 to b breaks
#' @param formula formula y~x object
#' @param data data
#' @param b number of breaks
#' @param ssr ssr table
#' @export
#' @importFrom magrittr %$%
bic <- function(formula, data, b, ssr) {

  z <- model.matrix(formula, data = data)
  q <- ncol(z)

  ssrval <- unlist(purrr::map(0:b, ~min_ssr(.x, ssr) %$% global_ssr))

  bicval <- unlist(purrr::map(1:(b+1), function(x) {

    log(ssrval[x] / nrow(data)) + log(nrow(data)) * (x-1) * (q+1) / nrow(data)

  }))

  minb <- which.min(bicval)

  list(
    breaks = minb-1,
    bic = setNames(bicval, 0:b),
    ssr = setNames(ssrval, 0:b)
  )

}


#' Calculate the modified Schwarz criterion for 0 to b breaks
#' @param formula formula y~x object
#' @param data data
#' @param b number of breaks
#' @param ssr ssr table
#' @export
#' @importFrom magrittr %$%
lwz <- function(formula, data, b, ssr) {

  z <- model.matrix(formula, data = data)
  q <- ncol(z)

  ssrval <- unlist(purrr::map(0:b, ~min_ssr(.x, ssr) %$% global_ssr))

  lwzval <- unlist(purrr::map(1:(b+1), function(x) {

    log(ssrval[x] / (nrow(data) - x * q - x + 1 )) + ((x - 1) * (q + 1) * 0.299 * (log(nrow(data)))^(2.1)) / nrow(data)

  }))

  minl <- which.min(lwzval)

  list(
    breaks = minl-1,
    lwz = setNames(lwzval, 0:b),
    ssr = setNames(ssrval, 0:b)
  )

}
