
#' Compute convex scagnostic measures
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
#'   a_s <- scree(anscombe$x1, anscombe$y1)
#'   sc_convex(a_s)
#'   a_s <- scree(anscombe$x2, anscombe$y2)
#'   sc_convex(a_s)
#'   a_s <- scree(anscombe$x3, anscombe$y3)
#'   sc_convex(a_s)
#'   a_s <- scree(anscombe$x4, anscombe$y4)
#'   sc_convex(a_s)
#' @export
sc_convex <- function(x, y) UseMethod("sc_convex")

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
sc_convex.default <- function(x, y){
  sc <- scree(x, y)
  sc_convex.scree(sc)
}

#' @rdname sc_convex
#' @export
sc_convex.list <- function(chull, ahull){
  chull_area <- splancs::areapl(cbind(chull$x, chull$y))
  if (ahull$length > 0)
    ahull_area <- alphahull::areaahull(ahull)
  else
    ahull_area <- 0
  ahull_area / chull_area
}

#' Compute convex scagnostic measures
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   anscombe_tidy <- anscombe %>%
#'   pivot_longer(cols = everything(),
#'     names_to = c(".value", "set"),
#'     names_pattern = "(.)(.)")
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   a_s <- scree(anscombe$x1, anscombe$y1)
#'   sc_skinny(a_s)
#'   a_s <- scree(anscombe$x2, anscombe$y2)
#'   sc_skinny(a_s)
#'   a_s <- scree(anscombe$x3, anscombe$y3)
#'   sc_skinny(a_s)
#'   a_s <- scree(anscombe$x4, anscombe$y4)
#'   sc_skinny(a_s)
#' @export
sc_skinny <- function(x, y) UseMethod("sc_skinny")

#' @rdname sc_skinny
#' @export
sc_skinny.scree <- function(x, y = NULL) {
  stopifnot(is.null(y))
  ahull <- gen_alpha_hull(x$del, x$alpha)
  sc_skinny.list(ahull)
}

#' @rdname sc_skinny
#' @export
sc_skinny.default <- function(x, y){
  sc <- scree(x, y)
  sc_skinny.scree(sc)
}

#' @rdname sc_skinny
#' @export
sc_skinny.list <- function(ahull){
  if (ahull$length > 0) {
    ahull_area <- alphahull::areaahull(ahull)
    s <- 1 - sqrt(4*pi * ahull_area) / ahull$length
  }
  else
    s <- 1
  return(s)
}

gen_conv_hull <- function(del) {
  tripack::convex.hull(del$tri.obj)
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
