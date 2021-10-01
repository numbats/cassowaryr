#' Measure of Discreteness
#'
#' This metric computes the ration between the
#' number of unique values to total data values
#' on number of rotations of the data, and
#' returns the smallest value. If this value is
#' small it means that there are only a few unique
#' data values, and hence the distribution is
#' discrete
#'
#' @param x numeric vector
#' @param y numeric vector
#' @return double
#' @examples
#' data("datasaurus_dozen_wide")
#' sc_striped(datasaurus_dozen_wide$v_lines_x,
#'            datasaurus_dozen_wide$v_lines_y)
#' sc_striped(datasaurus_dozen_wide$dino_x,
#'            datasaurus_dozen_wide$dino_y)
#' @export
sc_striped <- function(x, y){
  # TODO: Need to make this work for rotations,
  # and do some rounding
  dx <- unique(x)
  value <- 1-length(dx)/length(x)
  xs <- x - mean(x)
  ys <- y - mean(y)
  for (ang in seq(0, pi, pi/180)) {
    xr <- sin(ang)*xs + cos(ang)*ys
    #yr <- -cos(ang)*xs + sin(ang)*ys
    dx <- unique(xr)
    new_value <- 1-length(dx)/length(x)
    if (new_value > value) {
      value <- new_value
      cat(value, " ", new_value, "\n")
    }
  }
  value
}

