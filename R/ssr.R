#' Function to decide on valid segments in series
#' @param start begin of segment
#' @param end end of segment
#' @param length length of total observation period
#' @param n_breaks number of breaks in data
#' @param min_length minimum length of segment
#' @export
valid_segm <- function(start, end, length, n_breaks, min_length) {

  start <= end &
    ((start-1) %/% min_length + (length-end) %/% min_length) >= n_breaks &
    (end - start + 1) >= min_length &
    ((start-1) %/% min_length > 0 | start == 1)

}


#' calculate the squared recurive residual in a segment
#' @param formula y~x formula object
#' @param start start of segment
#' @param end end of segment
#' @param data dataframe of observations
#' @export
rec_res <- function(formula, start, end, data) {

  model <- lm(formula, data = data[start:end,])

  list(
    end = start:end,
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

  ssr_mtx <- matrix(nrow = nrow(cci), ncol = nrow(cci))

  min_length <- if(trim<1) floor(trim * nrow(data)) else max(trim, floor(0.05 * nrow(data)))

  df <- tibble::tibble(
    start = c(row(ssr_mtx)),
    end = c(col(ssr_mtx)),
    ssr = mapply(valid_segm, row(ssr_mtx), col(ssr_mtx), nrow(data), 0, min_length)
  )

  # calculate residuals necessary for sum squared residuals, avoiding recursive calculation
  dt <- df %>%
    dplyr::filter(ssr) %>%
    dplyr::group_by(start) %>%
    dplyr::slice_max(end) %>%
    dplyr::ungroup() %>%
    data.table::as.data.table()

  dt <- dt[,
           rec_res(formula, start, end, data),
           by = c("start")
  ]

  # SSR for start to end sequence
  residuals <- dt %>%
    tibble::as_tibble() %>%
    dplyr::group_by(start) %>%
    dplyr::mutate(
      ssr_val = cumsum(res_sq)
    ) %>%
    dplyr::ungroup()

  # join SSR with valid segments
  ssr_df <- df %>%
    dplyr::filter(ssr) %>%
    dplyr::left_join(
      residuals %>%
        dplyr::select(start, end, ssr_val),
      by = c("start", "end")
    )

  # create permutations of possible segment combinations
  ssr_df <- ssr_df %>%
    dplyr::mutate(
      end_start = start-1
    )

  return(ssr_df %>%
    dplyr::select(end, end_start, ssr = ssr_val)
  )

}


#' Dynamic programming algorithm for the global minimum of SSR for a number of breaks
#' @param breaks number of breaks
#' @param ssr ssr table
#' @export
#' @importFrom magrittr %>%
min_ssr <- function(breaks, ssr) {

  # dataframe of permutations
  permutations <- tibble::tibble(
    end_0 = 0
  )

  # join everything together at the breakpoints
  # and filter after one-break for minimal global ssr before next iteration
  # same as recursive solution but not recursive
  for(i in 1:(breaks+1)) {

    permutations <- permutations %>%
      dplyr::left_join(ssr %>%
                         dplyr::rename(setNames(c("end", "ssr"), c(paste0("end_", i), paste0("ssr_", i)))),
                       by = setNames("end_start", paste0("end_", i-1))) %>%
      dplyr::mutate(global_ssr = psum(!!!rlang::syms(paste("ssr", 1:i, sep = "_"))))

    if(i > 2) {
      permutations <- permutations %>%
        dplyr::group_by(!!rlang::sym(paste("end", i, sep = "_"))) %>%
        dplyr::slice_min(global_ssr, n = 1, with_ties = T) %>%
        dplyr::ungroup()
    }

  }

  # filter only permutation across full observation period
  permutations <- permutations %>%
    dplyr::slice_max(!!rlang::sym(paste0("end_", breaks+1))) %>%
    dplyr::slice_min(global_ssr)

  segments <- dplyr::bind_cols(
    permutations %>%
      dplyr::select(dplyr::all_of(paste("end", seq_len(breaks + 1), sep = "_"))) %>%
      tidyr::pivot_longer(cols = everything(), values_to = "end") %>%
      dplyr::select(end),
    permutations %>%
      dplyr::select(dplyr::all_of(paste("ssr", seq_len(breaks + 1), sep = "_"))) %>%
      tidyr::pivot_longer(cols = everything(), values_to = "ssr") %>%
      dplyr::select(ssr)
  ) %>%
    dplyr::mutate(start = dplyr::lag(end + 1, default = 1),
                  id = dplyr::row_number())

  list(
    breaks = breaks,
    start = segments$start,
    end = segments$end,
    id = segments$id,
    ssr = sum(segments$ssr)
  )
}
