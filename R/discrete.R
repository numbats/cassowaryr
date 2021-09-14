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
#' @export
sc_striped <- function(x, y){
  # TODO: Need to make this work for rotations,
  # and do some rounding
  dx <- unique(x)
  length(dx)/length(x)
}

