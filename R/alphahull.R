
#' Compute convex scagnostic measure
#'
#' @param x numeric vector of x values
#' @param y numeric vector of y values
#' @return A "numeric" object that gives the plot's convex score.
#'
#' @examples
#'   require(ggplot2)
#'   require(dplyr)
#'   ggplot(features, aes(x=x, y=y)) +
#'      geom_point() +
#'      facet_wrap(~feature, ncol = 5, scales = "free")
#'   features %>% group_by(feature) %>% summarise(convex = sc_convex(x,y))
#'   sc_convex(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#' @export
sc_convex <- function(x, y) UseMethod("sc_convex")

#' @rdname sc_convex
#' @export
sc_convex.default <- function(x, y){
  sc <- scree(x, y)
  sc_convex.scree(sc)
}

#' @rdname sc_convex
#' @export
sc_convex.scree <- function(x,y = NULL) {
  stopifnot(is.null(y))
  chull <- gen_conv_hull(x$del)
  ahull <- gen_alpha_hull(x$del, x$alpha)
  sc_convex.list(chull, ahull)
}

#' @rdname sc_convex
#' @export
sc_convex.list <- function(x, y){
  chull_area <- splancs::areapl(cbind(x$x, x$y))
  if (y$length > 0)
    ahull_area <- alphahull::areaahull(y)
  else
    ahull_area <- 0

  #calculate sample size weight
  #n = length(x)/500
  #w = 0.7 + 0.3/(1+n^2)
  (ahull_area / chull_area) # w*(ahull_area / chull_area)
}

#' Compute skinny scagnostic measure
#'
#' @param x numeric vector of x values
#' @param y numeric vector of y values
#' @return A "numeric" object that gives the plot's skinny score.
#'
#' @examples
#'   require(ggplot2)
#'   require(dplyr)
#'   ggplot(features, aes(x=x, y=y)) +
#'      geom_point() +
#'      facet_wrap(~feature, ncol = 5, scales = "free")
#'   features %>% group_by(feature) %>% summarise(skinny = sc_skinny(x,y))
#'   sc_skinny(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#' @export
sc_skinny <- function(x, y) UseMethod("sc_skinny")

#' @rdname sc_skinny
#' @export
sc_skinny.default <- function(x, y){
  sc <- scree(x, y)
  sc_skinny.scree(sc)
}

#' @rdname sc_skinny
#' @export
sc_skinny.scree <- function(x, y = NULL) {
  stopifnot(is.null(y))
  ahull <- gen_alpha_hull(x$del, x$alpha)
  sc_skinny.list(ahull)
}

#' @rdname sc_skinny
#' @export
sc_skinny.list <- function(x, y=NULL){
  if (x$length > 0) {
    ahull_area <- alphahull::areaahull(x)
    s <- 1 - sqrt(4*pi * ahull_area) / x$length
  }
  else
    s <- 1
  return(s)
}

gen_conv_hull <- function(del) {
  interp::convex.hull(del$tri.obj)
}

gen_alpha_hull <- function(del, alpha) {
  # This catches errors in ahull calculation and
  # returns a NULL if necessary

  ahull <- tryCatch(alphahull::ahull(del, alpha = alpha),
           error = function(c) {
             return(list(arcs=NULL, xahull=NULL, length=0, alpha=NULL))
           })
  return(ahull)
}
