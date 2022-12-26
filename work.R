
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
data <- read.table("tests/test.txt", header = T) %>% tibble::as_tibble()
data <- cci %>% dplyr::select(time_period, value)
data <- na.omit(cci %>% dplyr::select(time_period, diff))
data <- na.omit(cci %>% dplyr::select(time_period, diff_2))

# regression formula
fml <- as.formula(value ~ 1)
fml <- as.formula(diff ~ 1)
fml <- as.formula(diff_2 ~ 1)

# minimum number of observations in segment
trim <- 0.10

# build DF of SSR for later evaluation of different number of breaks
ssr <- ssr(fml, data, trim = trim)


# F-Statistics to find number of breaks

fstat <- purrr::map(1:4, ~fstats(fml, data, .x, ssr, trim = trim, autocorrelation = F))

dmax <- dmax_stats(fml, data, 4, ssr, trim = trim, autocorrelation = F)

lstat <- purrr::map(0:4, ~lstats(fml, data, .x, ssr, trim = trim, autocorrelation = F))

sequential <- seq_test(fml, data, 4, ssr, trim = trim, skip = -1, autocorrelation = F)

ssr_m <- purrr::map(0:4, ~min_ssr(.x, ssr))

bic <- bic(fml, data, 4, ssr)

lwz <- lwz(fml, data, 4, ssr)

scbp <- strucchange::breakpoints(data$diff_2 ~ 1, h = 0.1)
plot(scbp)
scbp_fs <- strucchange::Fstats(data$diff_2 ~ 1, from = 0.1)
plot(scbp_fs)
strucchange::sctest(scbp_fs)
#TODO BIC and LWZ

# global minimal SSR

ssr_result <- min_ssr(1, ssr)
