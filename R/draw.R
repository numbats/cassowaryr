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

#' Drawing the MST
#'
#' This function will draw the MST for a
#' scatterplot.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @examples
#' require(dplyr)
#' require(ggplot2)
#' data("features")
#' nl <- features %>% filter(feature == "nonlinear2")
#' draw_MST(nl$x, nl$y)
#' @export
draw_MST <- function(x, y, alpha=0.5) {
  scree <- scree(x, y)
  MST <- gen_mst(scree$del, scree$weights)
  xystartend <- as_tibble(scree[["del"]][["mesh"]])
  MST_mat <- twomstmat(MST, scree)$mat
  d_MST <- xystartend %>%
    group_by(ind1,ind2) %>%
    mutate(connected = ifelse(any(!MST_mat[ind1,ind2]==0), 1, 0)) %>%
    ungroup() %>%
    filter(connected==1)
  ggplot(d_MST) +
    geom_point(aes(x=x1, y=y1)) +
    geom_point(aes(x=x2, y=y2)) +
    geom_segment(aes(x=x1, xend=x2,
                     y=y1, yend=y2))
}

#' Drawing the Convex Hull
#'
#' This function will draw the Convex Hull for a
#' scatterplot.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @examples
#' data("features")
#' nl <- features %>% filter(feature == "nonlinear2")
#' draw_convexhull(nl$x, nl$y)
#' @export
draw_convexhull <- function(x, y, alpha=0.5) {

  #make scree and convex hull
  scree <- scree(x, y)
  chull <- gen_conv_hull(scree$del)

  #make data of start and end points of hull
  data <- tibble(x = sc$del$mesh[,3],
                 y = sc$del$mesh[,4])
  ends <- tibble(x1 = chull$x,
                 y1 = chull$y,
                 x2 = chull$x[c(2:length(chull$x),1)],
                 y2 = chull$y[c(2:length(chull$y),1)])

  #plot
  ggplot(data, aes(x,y)) +
    geom_point() +
    geom_segment(data=ends, aes(x=x1, y=y1, xend=x2, yend=y2))

}
