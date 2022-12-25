
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

data <- na.omit(cci %>% dplyr::select(time_period, diff_2))

# regression formula
fml <- as.formula(diff_2 ~ 1)

# minimum number of observations in segment
trim <- 0.05

# build DF of SSR for later evaluation of different number of breaks
ssr <- ssr(fml, data, trim = trim)


# F-Statistics to find number of breaks

fstat <- purrr::map(1:2, ~fstats(fml, data, .x, ssr, trim = trim))

dmax <- dmax_stats(fml, data, 1, ssr, trim = trim)

lstat <- purrr::map(0:2, ~lstats(fml, data, .x, ssr, trim = trim))

sequential <- seq_test(fml, data, 2, ssr, trim = trim, skip = -1)

ssr_m <- purrr::map(0:2, ~min_ssr(.x, ssr))

#TODO BIC and LWZ

# global minimal SSR

ssr_result <- min_ssr(1, ssr)
