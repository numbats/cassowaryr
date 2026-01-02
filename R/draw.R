#' Drawing the graph objects
#'
#' These functions will draw the graph objects that are used to compute the
#' scagnostics. They are useful for debugging and seeing the impact of
#' parameter adjustments, alpha bining, or outlier removal.
#' You can draw the MST, convex hull, and alpha hull
#' with each respective draw_* function.
#'
#' @inheritParams scree
#' @param fill set to TRUE if you want the polygon filled
#'
#' @return A ggplot object that shows the respective graph object
#'
#' @examples
#' require(dplyr)
#' require(ggplot2)
#' require(alphahull)
#'
#' cl <- features %>% filter(feature == "clusters")
#'
#' # draw the alpha hull
#' draw_alphahull(cl$x, cl$y)
#'
#' # draw the MST
#' draw_mst(cl$x, cl$y)
#'
#' # draw the convex hull
#' draw_convexhull(cl$x, cl$y)
#'
#' # You can utilise these functions to see the impact of parameter changes
#' draw_alphahull(cl$x, cl$y, alpha = "omega")
#' @name draw_functions

#' @rdname draw_functions
#' @export
draw_alphahull <- function(x, y, out.rm = FALSE,
                     binner = NULL, alpha = "rahman",
                     fill = FALSE){
  UseMethod("draw_alphahull")
}


#' @export
draw_alphahull.default <- function(x, y, out.rm = FALSE,
                             binner = NULL, alpha = "rahman",
                             fill = FALSE){
  # x1 <- x2 <- y1 <- y2 <- NULL
  sc <- scree(x, y, out.rm, binner, alpha = alpha)
  draw_alphahull.scree(sc, fill = fill)
}

#' @export
draw_alphahull.scree <- function(x, y=NULL, out.rm = FALSE,
                                 binner = NULL, alpha = "rahman",
                                 fill = FALSE){
  ahull <- alphahull::ahull(x$del, alpha=x$alpha)
  d_ahull <- tibble::tibble(x1 = ahull[["xahull"]][,1][ahull[["arcs"]][,7]],
                            y1 = ahull[["xahull"]][,2][ahull[["arcs"]][,7]],
                            x2 = ahull[["xahull"]][,1][ahull[["arcs"]][,8]],
                            y2 = ahull[["xahull"]][,2][ahull[["arcs"]][,8]])
  #scaled tibble
  d <- tibble::tibble(x=ahull[["xahull"]][,1], y=ahull[["xahull"]][,2])

  #make plot
  if(!fill){
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(data=d, ggplot2::aes(x, y),
                          colour = "black")+
      ggplot2::geom_segment(data=d_ahull,
                            ggplot2::aes(x=x1, xend=x2, y=y1, yend=y2))
  }
  if (fill){ #draws points on top of fill
    p <- ggplot2::ggplot() +
      ggplot2::geom_polygon(data=d_ahull,
                            ggplot2::aes(x=x1, y=y1),
                            alpha = 0.5) +
      ggplot2::geom_point(data=d, ggplot2::aes(x, y),
                          colour = "black")+
      ggplot2::geom_segment(data=d_ahull,
                            ggplot2::aes(x=x1, xend=x2, y=y1, yend=y2)
                            )
  }
  p

}

#' @rdname draw_functions
#' @export
draw_mst <- function(x, y, out.rm = FALSE, binner = NULL){
  UseMethod("draw_mst")
}

#' @export
draw_mst.default <- function(x, y, out.rm = FALSE,
                             binner = NULL){
  # ind1 <- ind2 <- connected <- x1 <- x2 <- y1 <- y2 <- NULL
  sc <- scree(x, y, out.rm, binner)
  draw_mst.scree(sc)
}

#' @export
draw_mst.scree <- function(x, y=NULL, out.rm = FALSE,
                           binner = NULL){
  mst <- gen_mst(x$del, x$weights)
  draw_mst.igraph(mst, x)
}

#' @export
draw_mst.igraph <- function(x, y, out.rm = FALSE,
                            binner = NULL){
  xystartend <- tibble::as_tibble(y[["del"]][["mesh"]])

  # get edge matrix
  n <- length(igraph::E(x)$weight) + 1
  MST_mat <- matrix(x[], nrow = n)
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

#' @export
#' @rdname draw_functions
draw_convexhull <- function(x, y, out.rm = FALSE, binner = NULL,
                            fill = FALSE){
  UseMethod("draw_convexhull")
}

#' @export
draw_convexhull.default <- function(x, y, out.rm = FALSE,
                                    binner = NULL, fill = FALSE){
  # x1 <- x2 <- y1 <- y2 <- NULL
  sc <- scree(x, y, out.rm, binner)
  draw_convexhull.scree(sc, fill = fill)
}

#' @export
draw_convexhull.scree <- function(x, y=NULL, out.rm = FALSE,
                                  binner = NULL, fill = FALSE){
  sc <- x
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
                                              xend=x2, yend=y2))

  if (fill){ # re-draws points on top of fill
    p <- p + ggplot2::geom_polygon(data=d_ends, alpha=0.5,
                                   ggplot2::aes(x=x1, y=y1)) +
    ggplot2::geom_segment(data = d_ends,
                          ggplot2::aes(x=x1, y=y1,
                                       xend=x2, yend=y2)) +
    ggplot2::geom_point(data = d,
                        ggplot2::aes(x, y))
  }
  p
  }
