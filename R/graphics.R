
plot_theme <- function(x) {
  list(
    if(x == "date") ggplot2::scale_x_date(breaks = scales::breaks_pretty(10))
    else ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(5)),
    ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(5)),
    ggplot2::theme_minimal(),
    ggplot2::theme(
      axis.text.x= ggplot2::element_text(angle=60, hjust=1),
      #panel.background = element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey50"),
      axis.line = ggplot2::element_line(colour = "grey50"),
      aspect.ratio = 3/5
    )
  )
}


#' Create Plots for Data
#' @param cci Data as defined in the cci() function
#' @export
#' @importFrom ggplot2 +
plots <- function(cci = cci()) {

  p <- list()

  p$full_ts <- ggplot2::ggplot(cci, ggplot2::aes(x=time_period, y=value)) +
    ggplot2::geom_line( color="orange", linewidth = 1) +
    ggplot2::xlab("Index Year") +
    ggplot2::ylab("Index Value") +
    plot_theme("date")

  p$diff_ts <- ggplot2::ggplot(na.omit(cci), ggplot2::aes(x=time_period, y=diff)) +
    ggplot2::geom_line( color="orange", linewidth = 1) +
    ggplot2::xlab("Index Year") +
    ggplot2::ylab("Index Change") +
    plot_theme("date")

  p$log_ts <- ggplot2::ggplot(cci, ggplot2::aes(x=time_period, y=log_val)) +
    ggplot2::geom_line( color="orange", linewidth = 1) +
    ggplot2::xlab("Index Year") +
    ggplot2::ylab("Index Logarithm") +
    plot_theme("date")

  p$log_diff_ts <- ggplot2::ggplot(na.omit(cci), ggplot2::aes(x=time_period, y=log_diff)) +
    ggplot2::geom_line(color="orange", linewidth = 1) +
    ggplot2::xlab("Index Year") +
    ggplot2::ylab("Index Logarithm Change") +
    plot_theme("date")

  p$annual_ts <- ggplot2::ggplot(cci, ggplot2::aes(x=time_period, y=annual_growth)) +
    ggplot2::geom_line( color="orange", linewidth = 1) +
    ggplot2::xlab("Index Year") +
    ggplot2::ylab("Index Annual Growth") +
    plot_theme("date")

  p$dens_ts <- ggplot2::ggplot(cci, ggplot2::aes(x = value)) +
    ggplot2::geom_density(colour = "orange", linewidth = 1) +
    ggplot2::xlab("Index Value") +
    ggplot2::ylab("density") +
    plot_theme("value")

  return(p)

}


#' Plot breaks in data
#' @param data cci dataframe
#' @param segments segments dataframe with start, end and id
#' @export
#' @importFrom magrittr %>%
plot_breaks <- function(data, segments) {

  p <- data %>%
    dplyr::mutate(row_number = dplyr::row_number()) %>%
    dplyr::left_join(segments, by = character()) %>%
    dplyr::filter(start <= row_number, end >= row_number)

  ggplot2::ggplot(p, ggplot2::aes(x=time_period, y=value, color = id)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::xlab("Index Year") +
    ggplot2::ylab("Index Value") +
    plot_theme("date") +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linetype = "dotted")

}
