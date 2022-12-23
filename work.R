
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



# break analysis ----------------------------------------------------------

# regression formula
fml <- as.formula(value ~ 1)

# minimum number of observations in segment
trim <- 0.10

# build DF of SSR for later evaluation of different number of breaks
ssr <- ssr(fml, cci, trim = trim)


# F-Statistics to find number of breaks

fstat <- purrr::map(1:7, ~fstats(fml, cci, .x, ssr, trim = trim))

dmax <- dmax_stats(fml, cci, 7, ssr, trim = trim)

lstat <- purrr::map(1:7, ~lstats(fml, cci, .x, ssr, trim = trim))

sequential <- seq_test(fml, cci, 6, ssr, trim = trim)

ssr_m <- purrr::map(1:7, ~min_ssr(.x, ssr))

#TODO BIC and LWZ

# global minimal SSR

ssr_result <- min_ssr(6, ssr)
