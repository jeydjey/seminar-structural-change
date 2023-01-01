
delta_coef <- function(formula, data, from, to) {
  #get coefficient for a segment
  lm(formula, data = data[from:to,])$coeff
}


vcovs <- function(formula, data, from, to, autocorrelation = FALSE) {
  #variance covariance matrix

  z <- model.matrix(formula, data)
  x <- model.response(model.frame(formula, data = data))

  #checking data distributions across segments
  data_dist_eq <- all(unlist(purrr::map(seq_len(length(from)-1), ~ks.test(x[from[.x]:to[.x]], x[from[.x+1]:to[.x+1]])$p.value))>=0.05)
  #checking error distribution across segments
  error_dist_eq <- all(unlist(purrr::map(seq_len(length(from)-1), ~ks.test(residuals(lm(formula, data[from[.x]:to[.x],])),
                                                                           residuals(lm(formula, data[from[.x+1]:to[.x+1],])))$p.value))>=0.05)
  zb <- zbar(z, from[-1]-1)

  if(!autocorrelation & !data_dist_eq & error_dist_eq){
    #No serial correlation, different distributions for the data, identical distributions for the errors
    mdl <- lm(formula, data)
    s_hat <- sigma(mdl)^2
    return(s_hat*solve((1/nrow(data)*t(zb)%*%zb)))
  }

  if(!autocorrelation & !data_dist_eq & !error_dist_eq){
    #No serial correlation, different distribution for the data, different distributions for the errors
    s_hat <- diag(
      unlist(purrr::map2(from, to, function(.x,.y) {
        d <- data[.x:.y,]
        mdl <- lm(formula, d)
        return(1/nrow(d) * sum(residuals(mdl)^2))
      }))
    )
    return(s_hat %x% diag(1, dim(z)[2]) %*% solve(t(zb) %*% zb))
  }

  if(!autocorrelation & data_dist_eq & error_dist_eq){
    #No serial correlation, identical distribution for the data, identical variances for the errors
    mdl <- lm(formula, data)
    s_hat <- sigma(mdl)^2
    l_hat <- diag(c(diff(from), nrow(data)-sum(diff(from))))/nrow(data)
    return(s_hat*solve(l_hat %*% (t(zb)%*%zb)))
  }

  if(!autocorrelation & data_dist_eq & !error_dist_eq){
    #No serial correlation, identical distribution for the data, different distributions for the errors
    s_hat <- diag(
      unlist(purrr::map2(from, to, function(.x,.y) {
        d <- data[.x:.y,]
        mdl <- lm(formula, d)
        return(1/nrow(d) * sum(residuals(mdl)^2))
      }))
    )
    l_hat <- diag(c(diff(from), nrow(data)-sum(diff(from))))/nrow(data)
    return(s_hat %x% diag(1, dim(z)[2]) %*% solve(l_hat %x% (t(zb)%*%zb)))
  }

  if(autocorrelation & !data_dist_eq & !error_dist_eq) {
    #serial correlation, different distributions for data and errors
    deltat <- diag(c(diff(from), nrow(data)-sum(diff(from))))
    hac <- purrr::map2(from, to, ~sandwich::kernHAC(lm(formula, data[.x:.y,]), prewhite = (.y-.x)>3))
    mhac <- matrix(0, nrow = dim(z)[2]*length(from), ncol = dim(z)[2]*length(from))
    mi <- seq(1, dim(z)[2]*length(from), dim(z)[2])
    for(i in 1:length(from)) {
      mhac[mi[i]:(mi[i]+1), mi[i]:(mi[i]+1)] <- drop(hac[[i]])
    }

    return(mhac)
  }

  if(autocorrelation & error_dist_eq) {
    #serial correlation, same distribution for the errors
    hac <- sandwich::kernHAC(lm(formula, data))
    l_hat <- diag(c(diff(from), nrow(data)-sum(diff(from))))/nrow(data)
    return(l_hat %x% hac)
  }

  stop(glue::glue("vcov case for autocorrelation = {autocorrelation}, equal data distribution = {data_dist_eq} and equal error distribution = {error_dist_eq} is not defined"),
      call. = FALSE)

}

fcrit <- function(conf, trim, q, k) {
  #critical f value
  dplyr::pull(dplyr::filter(readRDS(system.file("extdata", "crit_supf.rda", package = "seminar.bp.strucchange", mustWork = T)),
                conf == !!conf, trim == !!trim, q == !!q, k == !!k), crit)

}

dcrit <- function(conf, trim, q, dmax = c("udmax", "wdmax")) {
  #critical WDmax and UDmax value
  dplyr::pull(dplyr::filter(readRDS(system.file("extdata", "crit_dmax.rda", package = "seminar.bp.strucchange", mustWork = T)),
                            conf == !!conf, trim == !!trim, q == !!q, dmax == !!dmax), crit)
}

lcrit <- function(conf, trim, q, l) {
  #critical l+1|l value
  dplyr::pull(dplyr::filter(readRDS(system.file("extdata", "crit_seq.rda", package = "seminar.bp.strucchange", mustWork = T)),
                            conf == !!conf, trim == !!trim, q == !!q, l == !!l), crit)

}

#' Calculate F-statistic value
#' @param formula y~x formula object
#' @param data cci dataframe
#' @param breaks number of breaks to test for
#' @param ssr SSR Dataframe
#' @param conf confidence interval
#' @param trim trimming 0.05 - 0.25
#' @param autocorrelation does the data contain autocorrelation?
#' @export
fstats <- function(formula, data, breaks, ssr, conf = 0.95, trim = 0.1, autocorrelation = FALSE) {

  z <- model.matrix(formula, data = data)
  q <- ncol(z)

  obs <- nrow(data)

  #global minimum result for breaks number of breaks
  ssr_result <- min_ssr(breaks, ssr, init = min(ssr$to_from))

  #if no break is requested
  if(ssr_result$breaks == 0) {
    return(
      list(
        fstat = 0,
        breaks = 0,
        ssr = ssr_result$global_ssr,
        crit = NA_real_,
        significant = FALSE
      )
    )
  }

  #calculate coefs for segments
  deltas <- unlist(purrr::map2(ssr_result$from, ssr_result$to, ~delta_coef(formula, data, .x, .y)))
  mdelta <- matrix(deltas)

  #R matrix
  mr <- matrix(0, breaks, breaks+1)
  for(i in 1:breaks) {
    mr[i,i] <- -1
    mr[i,i+1] <- 1
  }
  mr <- mr %x% diag(q)

  #Variance covariance matrix
  mvcov <- vcovs(formula, data, ssr_result$from, ssr_result$to, autocorrelation = autocorrelation)

  #F statistic
  f <- ((obs-(breaks+1)*q)/(breaks*q*obs)) %*% t(mdelta) %*% t(mr) %*% solve(mr %*% mvcov %*% t(mr)) %*% mr %*% mdelta

  #Critical value
  crit <- fcrit(conf = conf, trim = trim, q = q, k = breaks)


  list(
    fstat = f[1,1],
    breaks = breaks,
    breaks_from = ssr_result$from,
    breaks_to = ssr_result$to,
    ssr = ssr_result$ssr,
    global_ssr = ssr_result$global_ssr,
    crit = crit,
    significant = f[1,1] > crit
  )

}

#' Double maximum test for up to the number of breaks
#' @param formula y~x formula object
#' @param data cci dataframe
#' @param breaks number of 1:breaks to test for
#' @param ssr SSR Dataframe
#' @param conf confidence interval
#' @param trim trimming 0.05 - 0.25
#' @param autocorrelation does the data contain autocorrelation?
#' @export
#' @importFrom magrittr %>% %$%
dmax_stats <- function(formula, data, breaks, ssr, conf = 0.95, trim = 0.1, autocorrelation = FALSE) {

  z <- model.matrix(formula, data = data)
  q <- ncol(z)

  #calculate fstats for 1:breaks breaks
  fstat <- dplyr::bind_rows(purrr::map(1:breaks, ~fstats(formula, data, .x, ssr, autocorrelation = autocorrelation) %$%
                               tibble::tibble(fstat = fstat, crit = crit)))

  #UD max
  ud <- max(fstat$fstat)

  #WD max
  wd <- fstat %>%
    dplyr::mutate(fstat = fstat * dplyr::first(crit)/crit) %>%
    dplyr::pull(fstat) %>%
    max()

  #Crit values
  ud_crit <- dcrit(conf = conf, trim = trim, q = q, dmax = "udmax")
  wd_crit <- dcrit(conf = conf, trim = trim, q = q, dmax = "wdmax")

  list(
    ud = ud,
    ud_crit = ud_crit,
    ud_significant = ud > ud_crit,
    wd = wd,
    wd_crit = wd_crit,
    wd_significant = wd > wd_crit
  )

}

#' F-test for l+1|l breaks
#' @param formula y~x formula object
#' @param data cci dataframe
#' @param breaks number of 1:breaks to test for
#' @param ssr SSR Dataframe
#' @param conf confidence interval
#' @param trim trimming 0.05 - 0.25
#' @param autocorrelation does the data contain autocorrelation?
#' @export
#' @importFrom magrittr %>%
lstats <- function(formula, data, l, ssr, conf = 0.95, trim = 0.1, autocorrelation = FALSE) {

  z <- model.matrix(formula, data = data)
  q <- ncol(z)

  #get global minimum ssr result for l breaks
  ssr_result <- min_ssr(l, ssr)

  # tests each segment for an additional break
  fstat <- dplyr::bind_rows(
    purrr::map2(ssr_result$from, ssr_result$to, function(.x, .y) {

      ssr_2 <- ssr %>%
        dplyr::filter(to <= .y, to_from >= (.x-1)) %>%
        dplyr::mutate(to = to - (.x-1),
                      to_from = to_from - (.x-1))

      fstats(formula, data[.x:.y,], 1, ssr_2, conf, trim, autocorrelation = autocorrelation)


    })
  )

  #max f stat
  fstat <- max(fstat$fstat)

  #critical value
  crit <- lcrit(conf = conf, trim = trim, q = q, l = l)

  list(
    l = l,
    fstat = fstat,
    crit = crit,
    significant = fstat > crit
  )

}

#' Sequential testing for another break up to m
#' @param formula y~x formula object
#' @param data cci dataframe
#' @param breaks number of 1:breaks to test for
#' @param ssr SSR Dataframe
#' @param conf confidence interval
#' @param trim trimming 0.05 - 0.25
#' @param skip Should the test of significance for skip+1|skip be skipped, -1 for no skip of 1|0
#' @param autocorrelation does the data contain autocorrelation?
#' @export
#' @importFrom magrittr %>%
seq_test <- function(formula, data, m, ssr, conf = 0.95, trim = 0.1, skip = -1, autocorrelation = FALSE) {

  z <- model.matrix(formula, data = data)
  q <- ncol(z)

  #global minimum without a break
  result <- min_ssr(0, ssr)

  #looping up to m breaks
  for(i in 1:m) {

    #find max fstat to insert new segment
    fstat <- purrr::map(1:length(result$from), function(.x) {

        .y <- list(
          from = result$from[.x],
          to = result$to[.x],
          id = result$id[.x],
          ssr = result$ssr[.x]
        )

        ssr_2 <- ssr %>%
          dplyr::filter(to <= .y$to, to_from >= (.y$from-1)) %>%
          dplyr::mutate(to = to - (.y$from-1),
                        to_from = to_from - (.y$from-1))

        stats <- fstats(formula, data[.y$from:.y$to,], 1, ssr_2, conf, trim, autocorrelation = autocorrelation)

        list(
          id = .y$id,
          fstat = stats$fstat,
          breaks_from = stats$breaks_from + .y$from - 1,
          breaks_to = stats$breaks_to + .y$from - 1,
          ssr = stats$ssr
        )

      }) %>% dplyr::bind_rows() %>%
      dplyr::slice_max(fstat)

    #check if !significance and return result unless skipping
    if(fstat$fstat[1] < lcrit(conf, trim, q, i-1) && !(i-1 <= skip))
      return(result)

    #insert break into segment
    old_id <- fstat$id[1]
    pos <- which(result$id==old_id)

    result$breaks <- result$breaks+1
    result$from <- c(result$from[-pos], fstat$breaks_from)
    o <- order(result$from)
    result$to <- c(result$to[-pos], fstat$breaks_to)[o]
    result$ssr <- c(result$ssr[-pos], fstat$ssr)[o]
    result$from <- result$from[o]
    result$id <- seq_len(result$breaks+1)
    result$global_ssr <- sum(result$ssr)

  }

  return(result)


}


#' Calculate Confidence intervals for 2.5% and 5% significance
#' @param formula formula y~x object
#' @param data data
#' @param from segment position from
#' @param to segment position to
#' @param autocorrellation autocorrelation present?
#' @export
conf_int <- function(formula, data, from, to, autocorrelation = FALSE) {

  y <- model.response(model.frame(formula, data = data))
  z <- model.matrix(formula, data = data)
  q <- ncol(z)

  #break positions
  breakpos <- from[-1] - 1

  zb <- zbar(z, breakpos)

  delta <- ols(y, zb)
  residuals <- y - zb %*% delta

  rb <- zbar(residuals, breakpos)

  #checking data distribution across segments
  data_dist_eq <- all(unlist(purrr::map(seq_len(length(from)-1), ~ks.test(y[from[.x]:to[.x]], y[from[.x+1]:to[.x+1]])$p.value))>=0.05)
  #checking error distribution across segments
  error_dist_eq <- all(unlist(purrr::map(seq_len(length(from)-1), ~ks.test(residuals(lm(formula, data[from[.x]:to[.x],])),
                                                                           residuals(lm(formula, data[from[.x+1]:to[.x+1],])))$p.value))>=0.05)
  #checking regressor distribution across segments
  regressor_dist_eq <- all(unlist(
    lapply(1:q, function(a, z) {
      all(unlist(purrr::map(seq_len(length(from)-1), ~ks.test(z[from[.x]:to[.x],a],
                                                              z[from[.x+1]:to[.x+1],a])$p.value))>=0.05)
    }, z = z)
  ))

  #coeff matrix
  dhat_temp <- delta[(q+1):(q*length(from))] - delta[1:(q*length(from)-q)]
  dmat <- matrix(0,  ncol = q*length(breakpos), nrow = length(breakpos))
  for(i in 1:length(breakpos)) {
    dmat[i,i*q-1] <- dhat_temp[i*q-1]
    dmat[i,i*q] <- dhat_temp[i*q]
  }

  #T of segment
  deltat <- diag(c(diff(from), nrow(data)-sum(diff(from))))


  #Q Matrix
  if(regressor_dist_eq) {
    #regressor distribution equal
    qmat1 <- matrix(0, nrow = q*length(breakpos), ncol = q*length(breakpos))
    for(i in 1:length(breakpos)) {
      qmat1[(i*q-1):(i*q),(i*q-1):(i*q)] <- t(z) %*% z / dim(z)[1]
    }

    qmat2 <- qmat1

  }

  if(!regressor_dist_eq) {
    #different distribution of regressors
    z_temp <- t(zb[,1:(q*length(breakpos))]) %*% zb[,1:(q*length(breakpos))]
    for(i in 1:(length(from)-1)) {
      z_temp[(i*q-1):(i*q),(i*q-1):(i*q)] <- z_temp[(i*q-1):(i*q),(i*q-1):(i*q)] / deltat[i,i]
    }
    qmat1 <- z_temp

    z_temp <- t(zb[,(q+1):(q*(length(breakpos)+1))]) %*% zb[,(q+1):(q*(length(breakpos)+1))]
    for(i in 1:(length(from)-1)) {
      z_temp[(i*q-1):(i*q),(i*q-1):(i*q)] <- z_temp[(i*q-1):(i*q),(i*q-1):(i*q)] / deltat[i+1,i+1]
    }
    qmat2 <- z_temp

  }


  #Phi
  if(!autocorrelation & error_dist_eq) {
    # no autocorrealation and residuals follow same distribution
    phimat1 <- matrix(0, nrow = length(breakpos), ncol = length(breakpos))
    for(i in 1:length(breakpos)) {
      phimat1[i,i] <- t(residuals) %*% residuals / dim(z)[1]
    }
    phimat2 <- phimat1
  }

  if(!autocorrelation & !error_dist_eq) {
    # no autocorrealation and residuals follow different distribution
    r_temp <- t(rb[,1:length(breakpos)]) %*% rb[,1:length(breakpos)]
    for(i in 1:length(breakpos)) {
      r_temp[i,i] <- r_temp[i,i] / deltat[i,i]
    }
    phimat1 <- r_temp

    r_temp <- t(rb[,2:(length(breakpos)+1)]) %*% rb[,2:(length(breakpos)+1)]
    for(i in 1:length(breakpos)) {
      r_temp[i,i] <- r_temp[i,i] / deltat[i+1,i+1]
    }
    phimat2 <- r_temp
  }

  if(autocorrelation & !error_dist_eq) {
    # autocorrealtion and residuals follow different distribution
    # independent HAC estimator for each segment
    hac <- purrr::map2(from, to, ~sandwich::kernHAC(lm(formula, data[.x:.y,]), prewhite = (.y-.x)>3))
    omega <- matrix(0, nrow = dim(z)[2]*length(from), ncol = dim(z)[2]*length(from))
    mi <- seq(1, dim(z)[2]*length(from), dim(z)[2])
    for(i in 1:length(from)) {
      omega[mi[i]:(mi[i]+1), mi[i]:(mi[i]+1)] <- drop(hac[[i]])
    }

    omega1 <- omega[1:(q*length(breakpos)), 1:(q*length(breakpos))]
    omega2 <- omega[(1+q):(q*length(breakpos)+q), (1+q):(q*length(breakpos)+q)]

  }

  if(autocorrelation & error_dist_eq) {
    # autocorrealtion and residuals follow same distribution
    # one HAC estimator for all segments
    hac <- sandwich::kernHAC(lm(formula, data))
    omega1 <- matrix(0, nrow = q*length(breakpos), ncol = q*length(breakpos))
    for(i in 1:length(breakpos)) {
      omega1[(i*q-1):(i*q), (i*q-1):(i*q)] <- drop(hac)
    }

    omega2 <- omega1
  }

  if(autocorrelation) {

    phimat1 <- dmat %*% omega1 %*% t(dmat) / (dmat %*% qmat1 %*% t(dmat))
    phimat2 <- dmat %*% omega2 %*% t(dmat) / (dmat %*% qmat2 %*% t(dmat))

  }

  eta <- dmat %*% qmat2 %*% t(dmat) /  (dmat %*% qmat1 %*% t(dmat))

  #critical value calculation
  crit <- crit_val(diag(eta), diag(phimat1), diag(phimat2))

  if(!autocorrelation) {
    a = diag(dmat %*% qmat1 %*% t(dmat) %*% solve(phimat1))
  } else {
    a = diag((dmat %*% qmat1 %*% t(dmat)) %*% (dmat %*% qmat1 %*% t(dmat)) / (dmat %*% omega1 %*% t(dmat)))
  }

  #boundaries
  bound <- matrix(0, nrow = 4, ncol = length(breakpos))

  bound[1,] <- floor(breakpos - crit[4,]/a)
  bound[2,] <- ceiling(breakpos - crit[1,]/a)+1
  bound[3,] <- floor(breakpos - crit[3,]/a)
  bound[4,] <- ceiling(breakpos - crit[2,]/a)+1

  rownames(bound) <- c("2.5% lower", "2.5% upper", "5% lower", "5% upper")

  return(bound)

}

funcg <- function (x,bet,alph,b,deld,gam){
  #function copied from GAUSS program of Bai and Perron

  if(x <= 0){
    xb = bet * sqrt(abs(x))
    if (abs(xb) <= 30){
      g = -sqrt(-x/(2*pi)) * exp(x/8) - (bet/alph)*exp(-alph*x)*pnorm(-bet*sqrt(abs(x)))+
        ((2*bet*bet/alph)-2-x/2)*pnorm(-sqrt(abs(x))/2)
    }
    else{
      aa=log(bet/alph)-alph*x-xb^2/2-log(sqrt(2*pi))-log(xb)
      g=-sqrt(-x/(2*pi))*exp(x/8)-exp(aa)*pnorm(-sqrt(abs(x))/2)+
        ((2*bet*bet/alph)-2-x/2)*pnorm(-sqrt(abs(x))/2)
    }
  } else {
    xb=deld*sqrt(x)

    if (abs(xb) <= 30){

      g=1+(b/sqrt(2*pi))*sqrt(x)*exp(-b*b*x/8)+
        (b*deld/gam)*exp(gam*x)*pnorm(-deld*sqrt(x))+
        (2-b*b*x/2-2*deld*deld/gam)*pnorm(-b*sqrt(x)/2);
    }
    else{

      aa=log((b*deld/gam))+gam*x-xb^2/2-log(sqrt(2*pi))-log(xb);
      g=1+(b/sqrt(2*pi))*sqrt(x)*exp(-b*b*x/8)+
        exp(aa)+(2-b*b*x/2-2*deld*deld/gam)*pnorm(-b*sqrt(x)/2);
    }
  }
  return(g)
}


crit_val <- function(eta, phi1, phi2) {
  #function copied from GAUSS program of Bai and Perron

  a <- phi1/phi2
  gam <- ((phi2/phi1)+1)*eta/2
  b <- sqrt(phi1*eta/phi2)
  deld <- sqrt(phi2*eta/phi1)+b/2
  alph <- a * (1+a)/2
  bet <- (1+2*a)/2
  sig <- c(0.025, 0.05, 0.95, 0.975)

  crit <- purrr::pmap(
    list(gam, b, deld, alph, bet),
    function(gam, b, deld, alph, bet) {

      out <- double(4)

      for(i in 1:4) {
        upb = 2000
        lwb = -2000
        crit = 999999
        cct = 1

        while(abs(crit) >= 0.000001) {
          cct = cct + 1
          if (cct > 100){
            print('the procedure to get critical values for the break dates has reached the upper bound on the number of iterations. This may happens in the procedure cvg. The resulting confidence interval for this break date is incorrect')
            break
          }
          else{
            xx <- lwb + (upb-lwb)/2
            pval <- funcg(xx,bet,alph,b,deld,gam)
            crit <- pval - sig[i]
            if (crit <= 0) {
              lwb <- xx}
            else {
              upb <- xx}
          }
        }
        out[i] = xx
      }

      return(out)

    }
  )

 m <- do.call(cbind, crit)

 rownames(m) <- sig

 return(m)

}



