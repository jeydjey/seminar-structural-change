

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
  xts$log <- xts::xts(cci$log_val, cci$time_period)
  xts$log_diff <- na.omit(xts::xts(cci$log_diff, cci$time_period))
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
  k$diff_mean <- mean(cci$diff)
  k$diff_var <- var(cci$diff)
  k$diff_below0 <- length(cci$diff[which(cci$diff<0)])

  #Test for stationarity with Augmented Dickey-Fuller Test and Phillips-Perron Unit Root Test
  #observation values
  kpi$adf_value <- tseries::adf.test(cci$value)
  kpi$pp_value <- tseries::pp.test(cci$value,)
  #first difference
  kpi$adf_diff <- tseries::adf.test(na.omit(cci$diff))
  kpi$pp_diff <- tseries::pp.test(na.omit(cci$diff))
  #log values
  kpi$adf_log_value <- tseries::adf.test(cci$log_val)
  kpi$pp_log_value <- tseries::pp.test(cci$log_val)
  #log first difference
  kpi$adf_log_diff <- tseries::adf.test(na.omit(cci$log_diff))
  kpi$pp_log_diff <- tseries::pp.test(na.omit(cci$log_diff))

  #kurtosis
  kpi$kurtosis <- TSA::kurtosis(cci$value)

  #Skewness
  kpi$skewness <- TSA::skewness(cci$value)

  # autocorrelation
  # observation values
  k$acf_value(xts$value, lag.max = 10, plot = F)
  k$pacf_value(xts$value, lag.max = 10, plot = F)
  # first difference
  k$acf_diff(xts$diff, lag.max = 10, plot = F)
  k$pacf_diff(xts$diff, lag.max = 10, plot = F)
  # log values
  k$acf_log(xts$log, lag.max = 10, plot = F)
  k$pacf_log(xts$log, lag.max = 10, plot = F)
  #first difference log values
  k$acf_log_diff(xts$log_diff, lag.max = 10, plot = F)
  k$pacf_log_diff(xts$log_diff, lag.max = 10, plot = F)

  return(k)

}
