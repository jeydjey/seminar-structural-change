install.packages(c("magrittr", "dplyr", "httr", "ggplot2", "readr"))

library(magrittr)



#Download EU27 Construction Cost Index quarterly data from 2000 to 2022
raw_data <- httr::GET(
  url = "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/STS_COPI_Q/Q.CSTI+CSTL+CSTM+CSTO.F_CC11_X_CC113.NSA.I15+I10.EU27_2020/?format=SDMX-CSV&startPeriod=2000-Q1&endPeriod=2022-Q3&lang=en&label=both"
) %>% 
  httr::content(as = "text") %>% 
  readr::read_csv() %>% 
  dplyr::rename_with(~tolower(.x)) %>% 
  dplyr::select(unit, geo, time_period, obs_value, obs_flag)

#alternative copy of above data is available at
#raw_data <- tbd

#converting raw data

