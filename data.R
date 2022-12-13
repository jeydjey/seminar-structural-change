install.packages(c("magrittr", "dplyr", "httr", "ggplot2", "readr", "lubridate", "tidyr", "glue", "zoo", "tseries", "TSA", "dynlm",
                   "strucchange", "data.table"))

library(magrittr)
library(ggplot2)
library(data.table)

# ggplot prep

#' @param x x-axis scale type
plot_theme <- function(x) {
  list(
    if(x == "date") scale_x_date(breaks = scales::breaks_pretty(10))
    else scale_x_continuous(breaks = scales::breaks_pretty(5)),
    scale_y_continuous(breaks = scales::breaks_pretty(5)),
    theme_minimal(),
    theme(
      axis.text.x=element_text(angle=60, hjust=1),
      #panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = "grey50"),
      axis.line = element_line(colour = "grey50"), 
      aspect.ratio = 3/5
    )
  )
}



# data preparation --------------------------------------------------------

plots <- list()
kpi <- list()

#Download EU27 Construction Cost Index quarterly data from 2000 to 2022
raw_data <- httr::GET(
  url = "https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/STS_COPI_Q?format=JSON&lang=en&freq=Q&indic_bt=CSTI&cpa2_1=CPA_F41001_X_410014&s_adj=NSA&unit=I15&geo=EU27_2020&time=2000-Q1&time=2000-Q2&time=2000-Q3&time=2000-Q4&time=2001-Q1&time=2001-Q2&time=2001-Q3&time=2001-Q4&time=2002-Q1&time=2002-Q2&time=2002-Q3&time=2002-Q4&time=2003-Q1&time=2003-Q2&time=2003-Q3&time=2003-Q4&time=2004-Q1&time=2004-Q2&time=2004-Q3&time=2004-Q4&time=2005-Q1&time=2005-Q2&time=2005-Q3&time=2005-Q4&time=2006-Q1&time=2006-Q2&time=2006-Q3&time=2006-Q4&time=2007-Q1&time=2007-Q2&time=2007-Q3&time=2007-Q4&time=2008-Q1&time=2008-Q2&time=2008-Q3&time=2008-Q4&time=2009-Q1&time=2009-Q2&time=2009-Q3&time=2009-Q4&time=2010-Q1&time=2010-Q2&time=2010-Q3&time=2010-Q4&time=2011-Q1&time=2011-Q2&time=2011-Q3&time=2011-Q4&time=2012-Q1&time=2012-Q2&time=2012-Q3&time=2012-Q4&time=2013-Q1&time=2013-Q2&time=2013-Q3&time=2013-Q4&time=2014-Q1&time=2014-Q2&time=2014-Q3&time=2014-Q4&time=2015-Q1&time=2015-Q2&time=2015-Q3&time=2015-Q4&time=2016-Q1&time=2016-Q2&time=2016-Q3&time=2016-Q4&time=2017-Q1&time=2017-Q2&time=2017-Q3&time=2017-Q4&time=2018-Q1&time=2018-Q2&time=2018-Q3&time=2018-Q4&time=2019-Q1&time=2019-Q2&time=2019-Q3&time=2019-Q4&time=2020-Q1&time=2020-Q2&time=2020-Q3&time=2020-Q4&time=2021-Q1&time=2021-Q2&time=2021-Q3&time=2021-Q4&time=2022-Q1&time=2022-Q2"
) %>% 
  httr::content(as = "text", encoding = "UTF-8") %>% 
  jsonlite::fromJSON(T) %$%
  tibble::tibble(time_period = names(dimension$time$category$index), value = unlist(value))

#alternative copy of above data is available at
#raw_data <- httr::GET(
  # url = https://github.com/jeydjey/seminar-structural-change/blob/8e076881494003237338339e89356019ff2ada68/data/raw_eu27_cci_data.rds
  # ) %>% 
  # httr::content(as = "raw") %>% 
  # unserialize()

#converting raw data
# data <- raw_data %>% 
#   tidyr::separate(indic_bt, into = c("index", "index_name"), sep = ":", remove = T) %>% 
#   dplyr::transmute(
#     index,
#     index_name,
#     time_period = lubridate::yq(time_period),
#     value = obs_value
#   ) %>% 
#   dplyr::filter(time_period <= lubridate::as_date("2022-04-01")) %>% 
#   dplyr::group_by(index) %>% 
#   dplyr::arrange(time_period) %>% 
#   dplyr::mutate(
#     diff = value - dplyr::lag(value, n = 1),
#     pct = diff / dplyr::lag(value, n = 1),
#     log_val = log(value),
#     annual_growth = 400 * log(value / dplyr::lag(value))
#   ) %>% 
#   dplyr::ungroup() %>% 
#   tidyr::replace_na(list(diff = 0, pct = 0))

#construction cost index data
cci <- raw_data %>% 
  dplyr::mutate(
    time_period = lubridate::yq(time_period),
    diff = value - dplyr::lag(value, n = 1),
    pct = diff / dplyr::lag(value, n = 1),
    log_val = log(value),
    annual_growth = 400 * log(value / dplyr::lag(value))
  ) %>%
  tidyr::replace_na(list(diff = 0, pct = 0))
  
xts_value <- xts::xts(cci$value, cci$time_period)
xts_log <- xts::xts(cci$log_val, cci$time_period)
xts_ann_growth <- xts::xts(cci$annual_growth, cci$time_period)
cci_xts <- ts(cci$value, start = c(2000, 1), end = c(2022, 2), frequency = 4)
# summary statistics ---------------------------------
#plot for full time series
plots$full_ts <- ggplot(cci, aes(x=time_period, y=value)) +
  geom_line( color="orange", linewidth = 1) + 
  xlab("Index Year") +
  ylab("Index Value") +
  plot_theme("date")

plots$full_ts
#plot for diff
plots$diff_ts <- ggplot(cci, aes(x=time_period, y=diff)) +
  geom_line( color="orange", linewidth = 1) + 
  xlab("Index Year") +
  ylab("Index Change") +
  plot_theme("date")

plots$diff_ts
#plot for logarithm
plots$log_ts <- ggplot(cci, aes(x=time_period, y=log_val)) +
  geom_line( color="orange", linewidth = 1) + 
  xlab("Index Year") +
  ylab("Index Logarithm") +
  plot_theme("date")

plots$log_ts

#plot for annual growth
plots$annual_ts <- ggplot(cci, aes(x=time_period, y=annual_growth)) +
  geom_line( color="orange", linewidth = 1) + 
  xlab("Index Year") +
  ylab("Index Annual Growth") +
  plot_theme("date")

plots$annual_ts
#calculating summary statistics

kpi$num_observations <- nrow(cci)
kpi$overall_first <- dplyr::first(cci$value)
kpi$overall_last <- dplyr::last(cci$value)
kpi$overall_min <- min(cci$value)
kpi$overall_max <- max(cci$value)
kpi$overall_mean <- mean(cci$value)
kpi$overall_var <- var(cci$value)
kpi$overall_change <- dplyr::last(cci$value) / dplyr::first(cci$value)
kpi$diff_mean <- mean(cci$diff)
kpi$diff_var <- var(cci$diff)
kpi$diff_below0 <- cci$diff[which(cci$diff<0)] %>% length()

#Augmented Dickey-Fuller Test
kpi$adf <- tseries::adf.test(cci$value)

#kurtosis
kpi$kurtosis <- TSA::kurtosis(cci$value)

#Skewness
kpi$skewness <- TSA::skewness(cci$value)

#plot for histogram

plots$dens_ts <- ggplot(cci, aes(x = value)) +
  geom_density(colour = "orange", linewidth = 1) +
  xlab("Index Value") +
  ylab("density") +
  plot_theme("value")

decompose(cci_xts, type = "additive") %>% plot()

# autocorrelation

acf(na.omit(xts_ann_growth), lag.max = 4, plot = F)


# break calculation -------------------------------------------------------

#' Function to decide on valid segments in series
#' @param start begin of segment
#' @param end end of segment
#' @param length length of total observation period
#' @param n_breaks number of breaks in data
#' @param min_length minimum length of segment
valid_segm <- function(start, end, length, n_breaks, min_length) {
  
  ((start-1) %/% min_length + (length-end) %/% min_length) >= n_breaks & 
    (end - start + 1) >= min_length &
    ((start-1) %/% min_length > 0 | start == 1) #& 
    #((length-end) %/% min_length > 0 | end == length)
  
}


#' calculate the nth squared residual in 1:n sequence
#' @param start start of segment
#' @param end end of segment
#' @param data vector of observations
rec_res <- function(start, end, data) {
  
  if(length(start)>1)
    return(purrr::map2(start, end, ~rec_res(.x, .y, data)) %>% unlist())
  
  
  model <- lm(data[start:end] ~ 1)
  
  v <- broom::augment(model) %>% 
    dplyr::select(.resid) %>% 
    dplyr::slice_tail(n = 1) %>% 
    dplyr::pull()
  
  return(v^2)
  
}

rec_res <- function(start, end, data) {
  
  model <- lm(data[start:end] ~1)
  
  list(
    end = start:end,
    res_sq = c(0, strucchange::recresid(model)^2)
  )
}

#' calculate a parallel sum
#' @param ... vectors to summarise
#' @param na.rm 
psum <- function(..., na.rm = FALSE) {
  dots <- rlang::dots_list(...)
  
  m <- do.call(cbind, dots)
  
  rowSums(m, na.rm = na.rm)
  
}


#' build empty matrix
ssr_mtx <- matrix(nrow = nrow(cci), ncol = nrow(cci))

#' set variables
#' minimal length 
min_length <- 6
#' number of breaks
n_breaks <- 5

#' build dataframe based on matrix size and create valid segments
ssr_df <- tibble::tibble(
  start = c(row(ssr_mtx)),
  end = c(col(ssr_mtx)),
  ssr = mapply(valid_segm, row(ssr_mtx), col(ssr_mtx), nrow(cci), n_breaks, min_length)
)

#' calculate residuals necessary for sum squared residuals, avoiding recursive calculation
#' 
dt <- ssr_df %>% 
  dplyr::filter(ssr) %>% 
  dplyr::group_by(start) %>% 
  dplyr::slice_max(end) %>% 
  dplyr::ungroup() %>% 
  data.table::as.data.table()

dt <- dt[,
  rec_res(start, end, cci$value),
  by = c("start")
]
  #purrr::map(~dplyr::mutate(.x, end = list(expand.grid(start, start:end)$Var2)) %>% tidyr::unnest(end)) %>% 
  #purrr::map(~dplyr::mutate(.x, res_sq = list(rec_res(start, end, cci$value))) %>% 
  #dplyr::bind_rows()

#' SSR for start to end sequence
residuals <- dt %>% 
  tibble::as_tibble() %>% 
  dplyr::group_by(start) %>% 
  dplyr::mutate(
    ssr_val = cumsum(res_sq)
  ) %>% 
  dplyr::ungroup()

#' join SSR with valid segments
ssr_df <- ssr_df %>% 
  dplyr::filter(ssr) %>% 
  dplyr::left_join(
    residuals %>% 
      dplyr::select(start, end, ssr_val),
    by = c("start", "end")
  )

#' create permutations of possible segment combinations
ssr_df <- ssr_df %>% 
  dplyr::mutate(
    end_start = start-1
  )

ssr <- ssr_df %>% 
  dplyr::select(end, end_start, ssr = ssr_val)

#' dataframe of permutations
permutations <- tibble::tibble(
  end_0 = 0
)

#' join everything together at the breakpoints
#' and filter after one-break for minimal global ssr before next iteration
#' same as recursive solution but not recursive
for(i in 1:(n_breaks+1)) {
  
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

#' filter only permutation across full observation period
permutations <- permutations %>% 
  dplyr::filter(get("==")(!!rlang::sym(paste0("end_", n_breaks+1)), nrow(cci)))

segments <- dplyr::bind_cols(
  permutations %>% 
      dplyr::select(dplyr::all_of(paste("end", seq_len(n_breaks + 1), sep = "_"))) %>% 
      tidyr::pivot_longer(cols = everything(), values_to = "end") %>% 
      dplyr::select(end),
  permutations %>% 
      dplyr::select(dplyr::all_of(paste("ssr", seq_len(n_breaks + 1), sep = "_"))) %>% 
      tidyr::pivot_longer(cols = everything(), values_to = "ssr") %>% 
      dplyr::select(ssr)
  ) %>% 
  dplyr::mutate(start = dplyr::lag(end + 1, default = 1),
                id = paste("segment", dplyr::row_number(), sep = "_"))

plot_data <- cci %>% 
  dplyr::mutate(row_number = dplyr::row_number()) %>% 
  dplyr::left_join(segments, by = character()) %>% 
  dplyr::filter(start <= row_number, end >= row_number)

ggplot(plot_data, aes(x=time_period, y=value, color = id)) +
  geom_line(linewidth = 1) + 
  xlab("Index Year") +
  ylab("Index Value") +
  plot_theme("date") +
  geom_smooth(method = "lm", se = FALSE)


