#' Drawing the alphahull
#'
#' This function will draw the alphahull for a
#' scatterplot.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @examples
#' data("features")
#' nl <- features %>% filter(feature == "nonlinear2")
#' draw_alphahull(nl$x, nl$y)
#' @export
draw_alphahull <- function(x, y, alpha=0.5) {
  d_ahull <- ahull(x, y, a=alpha)
  p <- ggplot(as_tibble(x, y), aes(x, y)) + geom_point()

  d_ahull_c <- d_ahull$ashape.obj
  p + geom_segment(data=as_tibble(d_ahull_c$edges),
                   aes(x=x1, xend=x2, y=y1, yend=y2))
}


