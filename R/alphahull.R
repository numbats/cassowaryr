
#' Compute convex scagnostic measure
#'
#' A measure of how convex the shape of the data is. It was first defined in
#' Graph Theoretic Scagnostics, Wilkinson, et al. (2005).Computed as the ratio
#' between the area of the alpha hull and convex hull. Unlike the other
#' scagnostic measures, a high value on convex does not correlate to an
#' interesting scatter plot, rather it usually indicates a lack of relationship
#' between the two variables.
#'
#' @inheritParams scree
#' @return A numeric object that gives the plot's convex score.
#' @examples
#'   require(ggplot2)
#'   require(dplyr)
#'   ggplot(features, aes(x=x, y=y)) +
#'      geom_point() +
#'      facet_wrap(~feature, ncol = 5, scales = "free")
#'   features |> group_by(feature) |> summarise(convex = sc_convex(x,y))
#'   sc_convex(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#' @export
sc_convex <- function(x, y, alpha = "rahman",
                      out.rm = TRUE, binner =  "hex") UseMethod("sc_convex")


#' @export
sc_convex.default <- function(x, y, alpha = "rahman",
                              out.rm = TRUE, binner =  "hex"){
  sc <- scree(x, y, alpha = alpha, out.rm = out.rm, binner = binner)
  sc_convex.scree(sc)
}


#' @export
sc_convex.scree <- function(x,y = NULL, alpha = "rahman",
                            out.rm = FALSE, binner = NULL) {
  stopifnot(is.null(y))
  chull <- gen_conv_hull(x$del)
  ahull <- gen_alpha_hull(x$del, x$alpha)
  sc_convex.list(chull, ahull)
}


#' @export
sc_convex.list <- function(x, y, alpha = "rahman",
                           out.rm = TRUE, binner =  "hex"){
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
#' A measure of how “thin” the shape of the data is. It was first defined in
#' Graph Theoretic Scagnostics, Wilkinson, et al. (2005). It is calculated
#' as the ratio between the area and perimeter of the alpha hull with some
#' normalisation such that 0 correspond to a perfect circle and values close
#' to 1 indicate a skinny polygon.
#'
#' @inheritParams scree
#' @return A numeric object that gives the plot's skinny score.
#' @examples
#'   require(ggplot2)
#'   require(dplyr)
#'   ggplot(features, aes(x=x, y=y)) +
#'      geom_point() +
#'      facet_wrap(~feature, ncol = 5, scales = "free")
#'   features |> group_by(feature) |> summarise(skinny = sc_skinny(x,y))
#'   sc_skinny(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#' @export
sc_skinny <- function(x, y, alpha = "rahman",
                      out.rm = TRUE, binner =  "hex") UseMethod("sc_skinny")

#' @export
sc_skinny.default <- function(x, y, alpha = "rahman",
                              out.rm = TRUE, binner =  "hex"){
  sc <- scree(x, y, alpha = alpha, out.rm = out.rm, binner = binner)
  if (is.null(sc$del)) {
    dc <- sc_dcor(x,y)
    if (dc >= 1 - 1e-8) {
      return(1)
    }

  }
  sc_skinny.scree(sc)
}


#' @export
sc_skinny.scree <- function(x, y = NULL, alpha = "rahman",
                            out.rm = FALSE, binner = NULL) {
  stopifnot(is.null(y))
  ahull <- gen_alpha_hull(x$del, x$alpha)
  sc_skinny.list(ahull)
}


#' @export
sc_skinny.list <- function(x, y=NULL, alpha = "rahman",
                           out.rm = FALSE, binner = NULL){
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

