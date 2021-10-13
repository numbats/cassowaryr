#' Measure of Spearman Correlation
#'
#' @param x numeric vector
#' @param y numeric vector
#' @return double
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe)
#'   anscombe_tidy <- anscombe %>%
#'   pivot_longer(cols = everything(),
#'     names_to = c(".value", "set"),
#'     names_pattern = "(.)(.)")
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_monotonic(anscombe$x1, anscombe$y1)
#'   sc_monotonic(anscombe$x2, anscombe$y2)
#'   sc_monotonic(anscombe$x3, anscombe$y3)
#'   sc_monotonic(anscombe$x4, anscombe$y4)
#' @export
sc_monotonic <- function(x, y){
  abs(stats::cor(x, y, method='spearman'))
}

#' Spline based index.
#'
#' (Taken from tourr git repo)
#' Compares the variance in residuals of a fitted
#' spline model to the overall variance to find
#' functional dependence in 2D projections
#' of the data.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe)
#'   anscombe_tidy <- anscombe %>%
#'   pivot_longer(cols = everything(),
#'     names_to = c(".value", "set"),
#'     names_pattern = "(.)(.)")
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_splines(anscombe$x1, anscombe$y1)
#'   sc_splines(anscombe$x2, anscombe$y2)
#'   sc_splines(anscombe$x3, anscombe$y3)
#' @export
sc_splines <- function(x,y) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package: mgcv must be installed to use splines scagnostic")
  }
  nx <- length(unique(x[!is.na(x)]))
  ny <- length(unique(y[!is.na(y)]))
  if (nx < 4 || ny < 4) {
    measure <- 0
  }
  else {
    kx <- ifelse(nx < 20, 3, 10)
    ky <- ifelse(ny < 20, 3, 10)
    mgam1 <- mgcv::gam(y ~ s(x, bs = "cr", k = kx))
    mgam2 <- mgcv::gam(x ~ s(y, bs = "cr", k = ky))
    measure <- max(1 - stats::var(stats::residuals(mgam1), na.rm = T) / stats::var(y, na.rm = T), 1 - stats::var(stats::residuals(mgam2), na.rm = T) / var(x, na.rm = T))
  }
  return(measure)
}

#' Distance correlation index.
#'
#' (Taken from tourr git repo)
#' Computes the distance correlation based index on
#' 2D projections of the data.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe)
#'   anscombe_tidy <- anscombe %>%
#'   pivot_longer(cols = everything(),
#'     names_to = c(".value", "set"),
#'     names_pattern = "(.)(.)")
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_dcor(anscombe$x1, anscombe$y1)
#'   sc_dcor(anscombe$x2, anscombe$y2)
#'   sc_dcor(anscombe$x3, anscombe$y3)
#'   sc_dcor(anscombe$x4, anscombe$y4)
#' @export
sc_dcor <- function(x,y) {
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop("Package: energy is required to compute dcor scagnostic.")
  }
  xy <- stats::na.omit(data.frame(x = x, y = y))
  measure <- with(xy, energy::dcor(x, y))
  return(measure)
}
