#' Measure of Spearman Correlation
#'
#' The measure of Spearman correlation calculated using the stats package
#' cor function with method='spearman'.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @return A "numeric" object that gives the plot's monotonic score.
#'
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
#'
#' @export
#' @seealso stats::cor
sc_monotonic <- function(x, y){
  abs(stats::cor(x, y, method='spearman'))
}

#' Spline based index.
#'
#' Measures the functional non-linear dependence by fitting a penalised
#' splines model on X using Y, and on Y using X. The measure was
#' defined as an association scagnostic in Katrin Grimm's PhD thesis (2016)
#' The variance of the residuals are scaled down by the axis so they are
#' comparable, and finally the maximum is taken.
#' Therefore the value will be closer to 1 if either relationship
#' can be decently explained by a splines model.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @return A "numeric" object that gives the plot's spines score.
#'
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
#'
#' @export
sc_splines <- function(x,y) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package: mgcv must be installed to use splines scagnostic")
  }
  # Check for near unique small numbers
  bx <- cut(x, 50)
  by <- cut(y, 50)
  nx <- length(unique(bx))
  ny <- length(unique(by))
  #nx <- length(unique(x[!is.na(x)]))
  #ny <- length(unique(y[!is.na(y)]))
  if (nx < 10 || ny < 10) {
    measure <- 0
  }
  else {
    kx <- ifelse(nx < 30, 3, 10) # Might need to tweak these some more
    ky <- ifelse(ny < 30, 3, 10) # number of knots used should be small
    mgam1 <- mgcv::gam(y ~ s(x, bs = "cr", k = kx))
    mgam2 <- mgcv::gam(x ~ s(y, bs = "cr", k = ky))
    measure <- max(1 - stats::var(stats::residuals(mgam1),
                                  na.rm = T) / stats::var(y, na.rm = T),
                   1 - stats::var(stats::residuals(mgam2), na.rm = T) /
                     stats::var(x, na.rm = T))
  }
  return(measure)
}

#' Distance correlation index.
#'
#' The distance correlation between X and Y defined by Székely, et al.
#' in Measuring and testing dependence by correlation of distances.
#' The measure was suggested as an  association scagnostic in
#' Katrin Grimm's PhD thesis (2016). Distance correlation is a measure
#' of non-linear dependence which is 0 if and only if the two variables
#' are independent. It is computed using an ANOVA like calculation on the
#' pairwise distances between observations.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @return A "numeric" object that gives the plot's dcor score.
#'
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
#'
#' @export
sc_dcor <- function(x,y) {
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop("Package: energy is required to compute dcor scagnostic.")
  }
  xy <- stats::na.omit(data.frame(x = x, y = y))
  measure <- with(xy, energy::dcor(x, y))
  return(measure)
}
