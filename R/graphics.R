
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
plots <- function(cci = cci()) {

  p <- list()

  p$full_ts <- ggplot2::ggplot(cci, ggplot2::aes(x=time_period, y=value)) +
    ggplot2::geom_line( color="orange", linewidth = 1) +
    ggplot2::xlab("Index Year") +
    ggplot2::ylab("Index Value") +
    plot_theme("date")

  p$diff_ts <- ggplot2::ggplot(na.omit(dplyr::select(cci, time_period, diff)), ggplot2::aes(x=time_period, y=diff)) +
    ggplot2::geom_line( color="orange", linewidth = 1) +
    ggplot2::xlab("Index Year") +
    ggplot2::ylab("Index First Difference") +
    plot_theme("date")

  p$diff_2_ts <- ggplot2::ggplot(na.omit(dplyr::select(cci, time_period, diff_2)), ggplot2::aes(x=time_period, y=diff_2)) +
    ggplot2::geom_line( color="orange", linewidth = 1) +
    ggplot2::xlab("Index Year") +
    ggplot2::ylab("Index Second Difference") +
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
#' @param value column of values to plot
#' @export
#' @importFrom magrittr %>%
plot_breaks <- function(data, segments, value = rlang::sym("value")) {

  val <- rlang::enquo(value)

  p <- data %>%
    dplyr::mutate(row_number = dplyr::row_number()) %>%
    dplyr::left_join(segments, by = character()) %>%
    dplyr::filter(from <= row_number, to >= row_number) %>%
    dplyr::mutate(Segment = as.character(id))

  ggplot2::ggplot(p, ggplot2::aes(x=time_period, y=!!val, color = Segment)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::xlab("Index Year") +
    ggplot2::ylab("Index Value") +
    ggplot2::scale_color_brewer(palette = "Set1") +
    plot_theme("date") +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linetype = "dotted")

}

#' annotate the plot
#' @param label label
#' @param x,y pos of data
#' @param xpos,ypos pos of label
#' @export
plot_annotate <- function(label, x, y, xpos, ypos, hjust) {

  curve <- dplyr::case_when(
    hjust == "left" & ypos >= y ~ .3,
    hjust == "right" & ypos >= y ~ -.3,
    hjust == "left" & ypos < y ~ -.3,
    hjust == "right" & ypos < y ~ .3
  )

  list(
    ggplot2::annotate(
      geom = "curve", x = xpos, y = ypos, xend = x, yend = y,
      curvature = curve, arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm"))),
    ggplot2::annotate(geom = "text", x = xpos+ ifelse(hjust == "left", 0.1, -0.1), y = ypos, label = label, hjust = hjust)
  )
}

#' draw confidence interval
#' @param x lower and upper bound of x
#' @param y y position
#' @export
plot_confidence <- function(x, y) {

  list(
    ggplot2::geom_errorbar(data = NULL, ggplot2::aes(y = y, xmin = x[1], xmax = x[2]), colour = "red",
                           show.legend = FALSE, inherit.aes = FALSE)
  )

}
