
# setup -------------------------------------------------------------------

#GitHub Repo @ https://github.com/jeydjey/seminar-structural-change

devtools::install_github("jeydjey/seminar-structural-change")

library(magrittr)
library(seminar.bp.strucchange)

cci <- cci()


# data analysis -----------------------------------------------------------

#load kpi
kpi <- kpi(cci)

#load plots
plots <- plots(cci)

plots$acf_value <- plot(kpi$acf_value, main = "Construction Cost Index values")
plots$acf_diff <- plot(kpi$acf_diff, main = "Construction Cost Index first difference")
plots$acf_diff_2 <- plot(kpi$acf_diff_2, main = "Construction Cost Index second difference")


# break analysis ----------------------------------------------------------
data <- cci %>% dplyr::select(time_period, value, diff, diff_2) %>% na.omit()

# regression formula
fml <- as.formula(value ~ diff)

# minimum number of observations in segment
trim <- 0.10

# build DF of SSR for later evaluation of different number of breaks
ssr <- ssr(fml, data, trim = trim)

# F-Statistics to find number of breaks

fstat <- purrr::map(1:6, ~fstats(fml, data, .x, ssr, trim = trim, autocorrelation = T))

dmax <- dmax_stats(fml, data, 6, ssr, trim = trim, autocorrelation = T)

lstat <- purrr::map(0:6, ~lstats(fml, data, .x, ssr, trim = trim, autocorrelation = T))

sequential <- seq_test(fml, data, 6, ssr, trim = trim, skip = 0, autocorrelation = T)

ssr_m <- purrr::map(0:6, ~min_ssr(.x, ssr))

bic <- bic(fml, data, 6, ssr)

lwz <- lwz(fml, data, 6, ssr)

# global minimal SSR

ssr_result <- sequential

# K-S Test on data and errors

z <- model.matrix(fml, data)
x <- model.response(model.frame(fml, data = data))

#checking data distributions across segments
data_dist_eq <- all(unlist(purrr::map(seq_len(length(ssr_result$from)-1), ~ks.test(x[ssr_result$from[.x]:ssr_result$to[.x]], x[ssr_result$from[.x+1]:ssr_result$to[.x+1]])$p.value))>=0.05)
#checking error distribution across segments
error_dist_eq <- all(unlist(purrr::map(seq_len(length(ssr_result$from)-1), ~ks.test(residuals(lm(fml, data[ssr_result$from[.x]:ssr_result$to[.x],])),
                                                                         residuals(lm(fml, data[ssr_result$from[.x+1]:ssr_result$to[.x+1],])))$p.value))>=0.05)

# Tests in pretty
#F-stat
purrr::map(fstat, ~rlist::list.exclude(.x, !.name %in% c("fstat", "breaks", "significant"))) %>% dplyr::bind_rows()
#L-stat
purrr::map(lstat, ~rlist::list.exclude(.x, !.name %in% c("fstat", "l", "significant"))) %>% dplyr::bind_rows()
#

# Parameters with 3 breaks

delta <- purrr::map2(ssr_result$from, ssr_result$to, function(a,b) {
  y <- matrix(model.response(model.frame(fml, data = data[a:b,])))
  mm <- model.matrix(fml, data = data[a:b,])
  z <- matrix(mm, ncol = dim(mm)[2], byrow = F)

  ols(y, z)
}
)

#confidence intervals
conf <- conf_int(fml, data, ssr_result$from, ssr_result$to, autocorrelation = T)[3:4,]

break1_conf <- data$time_period[conf[,1]]
break2_conf <- data$time_period[conf[,2]]
break3_conf <- data$time_period[conf[,3]]


# Plot with the result of 3 breaks
label <- data %>%
  dplyr::filter(dplyr::row_number() %in% ssr_result$to, dplyr::row_number() != 88) %>%
  dplyr::mutate(
    breaks = glue::glue("break {dplyr::row_number()}: {lubridate::year(time_period)}-{lubridate::quarter(time_period)}")
  ) %>%
  split(1:nrow(.)) %>%
  lapply(as.list)

plot_breaks(data, dplyr::bind_rows(ssr_result), value = value) +
  plot_annotate(label = label[[1]]$breaks, x = label[[1]]$time_period, y = label[[1]]$value, xpos = lubridate::as_date("2006-01-01"), ypos = 98, hjust = "right") +
  plot_annotate(label = label[[2]]$breaks, x = label[[2]]$time_period, y = label[[2]]$value, xpos = lubridate::as_date("2014-01-01"), ypos = 108, hjust = "right") +
  plot_annotate(label = label[[3]]$breaks, x = label[[3]]$time_period, y = label[[3]]$value, xpos = lubridate::as_date("2016-01-01"), ypos = 116, hjust = "right") +
  plot_annotate(label = paste("potential", "break 4:", glue::glue("{lubridate::year(data[81,]$time_period)}-{lubridate::quarter(data[81,]$time_period)}"), sep = " "), x = data[81,]$time_period, y = data[81,]$value, xpos = lubridate::as_date("2018-01-01"), ypos = 93, hjust = "right") +
  plot_confidence(break1_conf, 70) +
  plot_confidence(break2_conf, 70) +
  plot_confidence(break3_conf, 70)




