#' Measure of Discreteness
#'
#' This metric computes the 1-(ratio between the
#' number of unique values to total data values)
#' on number of rotations of the data, and
#' returns the smallest value. If this value is
#' large it means that there are only a few unique
#' data values, and hence the distribution is
#' discrete
#'
#' @param x numeric vector
#' @param y numeric vector
#'
#' @return double
#' @examples
#' data("datasaurus_dozen_wide")
#' sc_striped(datasaurus_dozen_wide$v_lines_x,
#'            datasaurus_dozen_wide$v_lines_y)
#' sc_striped(datasaurus_dozen_wide$dino_x,
#'            datasaurus_dozen_wide$dino_y)
#' @export
sc_striped <- function(x, y){
  # Put both variables in range of -1 to 1
  xs <- 2 * ((x - min(x, na.rm = T))/
    (max(x, na.rm = T) - min(x, na.rm = T)) - 0.5)
  ys <- 2 * ((y - min(y, na.rm = T))/
               (max(y, na.rm = T) - min(y, na.rm = T)) - 0.5)
  # Bin to effectively round real values
  bx <- cut(xs, 50)
  # Only bins with observations
  dx <- unique(bx)
  # This is the index value
  value <- 1-length(dx)/length(levels(bx))
  # Check for a full range of 2D rotations
  for (ang in seq(0, pi, pi/180)) {
    xr <- sin(ang)*xs + cos(ang)*ys
    bx <- cut(xr, 50)
    dx <- unique(bx)
    new_value <- 1-length(dx)/length(levels(bx))
    if (new_value > value) {
      value <- new_value
    }
  }
  value
}

