#' Drawing the alphahull
#'
#' This function will draw the alphahull for a
#' scatterplot.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param alpha transparency value of points
#' @param clr optional colour of points and lines, default black
#' @param fill Fill the polygon
#' @return A "gg" object that draws the plot's alpha hull.
#' @examples
#' require(dplyr)
#' require(ggplot2)
#' require(alphahull)
#' data("features")
#' nl <- features %>% filter(feature == "nonlinear2")
#' draw_alphahull(nl$x, nl$y)
#' @export
draw_alphahull <- function(x, y, alpha=0.2, clr = "black", fill = FALSE) {
  x1 <- x2 <- y1 <- y2 <- NULL
  d_ahull <- alphahull::ahull(x, y, a=alpha)
  d <- tibble::tibble(x=x, y=y)
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data=d, ggplot2::aes(x, y),
                 colour = "black", alpha=0.5)

  d_ahull_c <- d_ahull$ashape.obj
  p <- p + ggplot2::geom_segment(data=tibble::as_tibble(d_ahull_c$edges),
                   ggplot2::aes(x=x1, xend=x2, y=y1, yend=y2),
                   colour = clr)
  p
}

#' Drawing the MST
#'
#' This function will draw the MST for a
#' scatterplot.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param alpha The alpha value used to build the graph object. Larger values allow points further apart to be connected.
#' @return A "gg" object that draws the plot's MST.
#' @examples
#' require(dplyr)
#' require(ggplot2)
#' data("features")
#' nl <- features %>% filter(feature == "nonlinear2")
#' draw_mst(nl$x, nl$y)
#' @export
draw_mst <- function(x, y, alpha=0.5) {
  ind1 <- ind2 <- connected <- x1 <- x2 <- y1 <- y2 <- NULL
  scree <- scree(x, y)
  MST <- gen_mst(scree$del, scree$weights)
  xystartend <- tibble::as_tibble(scree[["del"]][["mesh"]])
  MST_mat <- twomstmat(MST, scree)$mat
  d_MST <- xystartend %>%
    dplyr::group_by(ind1,ind2) %>%
    dplyr::mutate(connected = ifelse(any(!MST_mat[ind1,ind2]==0), 1, 0)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(connected==1)
  ggplot2::ggplot(d_MST) +
    ggplot2::geom_point(ggplot2::aes(x=x1, y=y1), alpha=0.5) +
    ggplot2::geom_point(ggplot2::aes(x=x2, y=y2), alpha=0.5) +
    ggplot2::geom_segment(ggplot2::aes(x=x1, xend=x2,
                     y=y1, yend=y2))
}

#' Drawing the Convex Hull
#'
#' This function will draw the Convex Hull for a
#' scatterplot.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param alpha transparency value of points
#' @param clr optional colour of points and lines, default black
#' @param fill Fill the polygon
#' @return A "gg" object that draws the plot's convex hull.
#' @examples
#' require(dplyr)
#' require(ggplot2)
#' data("features")
#' nl <- features %>% filter(feature == "nonlinear2")
#' draw_convexhull(nl$x, nl$y)
#' @export
draw_convexhull <- function(x, y, alpha=0.5, clr = "black", fill = FALSE) {
  x1 <- x2 <- y1 <- y2 <- NULL # happy cran checks
  # make scree and convex hull
  sc <- scree(x, y)
  chull <- gen_conv_hull(sc$del)

  # make data of start and end points of hull
  d <- tibble::tibble(x = sc$del$x[,1],
                 y = sc$del$x[,2])
  d_ends <- tibble::tibble(x1 = chull$x,
                 y1 = chull$y,
                 x2 = chull$x[c(2:length(chull$x),1)],
                 y2 = chull$y[c(2:length(chull$y),1)])

  # plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = d, ggplot2::aes(x, y), colour = "black", alpha = 0.5)

  p <- p + ggplot2::geom_segment(data = d_ends,
              ggplot2::aes(x=x1, y=y1,
                     xend=x2, yend=y2), colour = clr)

  if (fill) # re-draws points on top of fill
    p <- p + ggplot2::geom_polygon(data=d_ends,
                                   ggplot2::aes(x=x1, y=y1),
                          fill = clr, alpha = 0.5) +
         ggplot2::geom_segment(data = d_ends,
                          ggplot2::aes(x=x1, y=y1,
                                       xend=x2, yend=y2),
                          colour = clr) +
         ggplot2::geom_point(data = d,
                             ggplot2::aes(x, y),
                             colour = "black", alpha = 0.5)

  p
}
