

#' Read Construction Cost Index
#' @param import Import data from source?
#' @export
#' @importFrom magrittr %>%
cci <- function(import = FALSE) {

  if(!import) {
    return(readRDS(system.file("extdata", "cci.rda", package = "seminar.bp.strucchange", mustWork = T)))
  }

  raw_data <- httr::GET(
    url = "https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/STS_COPI_Q?format=JSON&lang=en&freq=Q&indic_bt=CSTI&cpa2_1=CPA_F41001_X_410014&s_adj=NSA&unit=I15&geo=EU27_2020&time=2000-Q1&time=2000-Q2&time=2000-Q3&time=2000-Q4&time=2001-Q1&time=2001-Q2&time=2001-Q3&time=2001-Q4&time=2002-Q1&time=2002-Q2&time=2002-Q3&time=2002-Q4&time=2003-Q1&time=2003-Q2&time=2003-Q3&time=2003-Q4&time=2004-Q1&time=2004-Q2&time=2004-Q3&time=2004-Q4&time=2005-Q1&time=2005-Q2&time=2005-Q3&time=2005-Q4&time=2006-Q1&time=2006-Q2&time=2006-Q3&time=2006-Q4&time=2007-Q1&time=2007-Q2&time=2007-Q3&time=2007-Q4&time=2008-Q1&time=2008-Q2&time=2008-Q3&time=2008-Q4&time=2009-Q1&time=2009-Q2&time=2009-Q3&time=2009-Q4&time=2010-Q1&time=2010-Q2&time=2010-Q3&time=2010-Q4&time=2011-Q1&time=2011-Q2&time=2011-Q3&time=2011-Q4&time=2012-Q1&time=2012-Q2&time=2012-Q3&time=2012-Q4&time=2013-Q1&time=2013-Q2&time=2013-Q3&time=2013-Q4&time=2014-Q1&time=2014-Q2&time=2014-Q3&time=2014-Q4&time=2015-Q1&time=2015-Q2&time=2015-Q3&time=2015-Q4&time=2016-Q1&time=2016-Q2&time=2016-Q3&time=2016-Q4&time=2017-Q1&time=2017-Q2&time=2017-Q3&time=2017-Q4&time=2018-Q1&time=2018-Q2&time=2018-Q3&time=2018-Q4&time=2019-Q1&time=2019-Q2&time=2019-Q3&time=2019-Q4&time=2020-Q1&time=2020-Q2&time=2020-Q3&time=2020-Q4&time=2021-Q1&time=2021-Q2&time=2021-Q3&time=2021-Q4&time=2022-Q1&time=2022-Q2"
  ) %>%
    httr::content(as = "text", encoding = "UTF-8") %>%
    jsonlite::fromJSON(T) %$%
    tibble::tibble(time_period = names(dimension$time$category$index), value = unlist(value))

  raw_data %>%
    dplyr::mutate(
      time_period = lubridate::yq(time_period),
      diff = value - dplyr::lag(value, n = 1),
      diff_2 = diff - dplyr::lag(diff, n = 1),
      pct = diff / dplyr::lag(value, n = 1),
      log_val = log(value),
      log_diff = log_val - dplyr::lag(log_val, n = 1),
      annual_growth = 400 * log(value / dplyr::lag(value))
    )

}


#' Return KPI for Construction Cost Index data
#'
#' @param cci Data as defined in the cci() function
#' @export
kpi <- function(cci = cci()) {

  xts <- list()

  xts$value <- xts::xts(cci$value, cci$time_period)
  xts$diff <- na.omit(xts::xts(cci$diff, cci$time_period))
  xts$diff_2 <- na.omit(xts::xts(cci$diff_2, cci$time_period))
  xts$ann_growth <- xts::xts(cci$annual_growth, cci$time_period)

  k <- list()

  k$num_observations <- nrow(cci)
  k$overall_first <- dplyr::first(cci$value)
  k$overall_last <- dplyr::last(cci$value)
  k$overall_min <- min(cci$value)
  k$overall_max <- max(cci$value)
  k$overall_mean <- mean(cci$value)
  k$overall_var <- var(cci$value)
  k$overall_change <- dplyr::last(cci$value) / dplyr::first(cci$value)
  k$diff_mean <- mean(cci$diff, na.rm = TRUE)
  k$diff_var <- var(cci$diff, na.rm = TRUE)
  k$diff_below0 <- length(cci$diff[which(cci$diff<0)])

  #Test for stationarity with Augmented Dickey-Fuller Test and Phillips-Perron Unit Root Test
  #observation values
  k$adf_value <- tseries::adf.test(cci$value)
  k$pp_value <- tseries::pp.test(cci$value,)
  #first difference
  k$adf_diff <- tseries::adf.test(na.omit(cci$diff))
  k$pp_diff <- tseries::pp.test(na.omit(cci$diff))
  #second difference
  k$adf_diff_2 <- tseries::adf.test(na.omit(cci$diff_2))
  k$pp_diff_2 <- tseries::pp.test(na.omit(cci$diff_2))

  #kurtosis
  k$kurtosis <- TSA::kurtosis(cci$value)

  #Skewness
  k$skewness <- TSA::skewness(cci$value)

  # autocorrelation
  # observation values
  k$acf_value <- acf(xts$value, lag.max = 10, plot = F)
  k$pacf_value <- pacf(xts$value, lag.max = 10, plot = F)
  k$dwt_value <- car::durbinWatsonTest(lm(value~1, cci), max.lag = 4)
  # first difference
  k$acf_diff  <- acf(xts$diff, lag.max = 10, plot = F)
  k$pacf_diff <- pacf(xts$diff, lag.max = 10, plot = F)
  k$dwt_diff <- car::durbinWatsonTest(lm(na.omit(cci$diff)~1), max.lag=4)
  #second difference
  k$acf_diff_2  <- acf(xts$diff_2, lag.max = 10, plot = F)
  k$pacf_diff_2 <- pacf(xts$diff_2, lag.max = 10, plot = F)
  k$dwt_diff_2 <- car::durbinWatsonTest(lm(na.omit(cci$diff_2)~1), max.lag=4)


  #evaluating model fit
  k$mdl_1 <- list(
    formula = value ~ 1,
    kpi = broom::glance(lm(value ~ 1, data = na.omit(cci)))
  )

  k$mdl_2 <- list(
    formula = value ~ diff,
    kpi = broom::glance(lm(value ~ diff, data = na.omit(cci)))
  )

  k$mdl_3 <- list(
    formula = value ~ diff + diff_2,
    kpi = broom::glance(lm(value ~ diff + diff_2, data = na.omit(cci)))
  )


  return(k)

}
