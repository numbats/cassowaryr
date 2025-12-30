#' Drawing the alphahull
#'
#' This function will draw the alphahull for a
#' scatterplot.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param fill Fill the polygon
#' @param binner an optional function that bins the x and y vectors prior
#' to triangulation
#'  Can be:
#'   \itemize{
#'     \item `NULL`: no binning (use raw points)
#'     \item `"hex"` : hexagonal binning following the procedure in the
#'     graph-theoretic scagnostics paper (start 40x40, halve
#'              until <= 250 nonempty cells)
#'     \item a function: user-defined binner
#'     }
#' @param alpha character, numeric, or function. Controls the alpha radius.
#'   Valid character values are:
#'   \itemize{
#'     \item "rahman" (default): Rahman's MST-based middle-50% alpha
#'     \item "q90": 90th percentile of MST edge lengths
#'     \item "omega": graph-theoretic scagnostics alpha
#'   Alternatively:
#'     \item a numeric value giving a fixed alpha
#'     \item a function with no arguments that returns a single numeric alpha
#'     }
#' @param outlier_rm logical; if TRUE, iteratively trim large MST edges
#' @return A alphahull::ahull(del, alpha = alpha) "gg" object that draws the plot's alpha hull.
#' @examples
#' require(dplyr)
#' require(ggplot2)
#' require(alphahull)
#' data("features")
#' cl <- features %>% filter(feature == "clusters")
#' draw_alphahull(cl$x, cl$y)
#' draw_alphahull(cl$x, cl$y, alpha = "omega")
#' @export
draw_alphahull <- function(x, y, outlier_rm = FALSE,
                     binner = NULL, alpha = "rahman",
                     fill = FALSE){
  UseMethod("draw_alphahull")
}

#' @rdname draw_alphahull
#' @export
draw_alphahull.default <- function(x, y, outlier_rm = FALSE,
                             binner = NULL, alpha = "rahman",
                             fill = FALSE){
  # x1 <- x2 <- y1 <- y2 <- NULL
  sc <- scree(x, y, outlier_rm, binner, alpha = alpha)
  draw_alphahull.scree(sc, fill = fill)
}

#' @rdname draw_alphahull
#' @export
draw_alphahull.scree <- function(x, y=NULL, outlier_rm = FALSE,
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

#' @export
draw_alphahull_old <- function(x, y, fill = FALSE, out.rm=TRUE) {
  x1 <- x2 <- y1 <- y2 <- NULL
  #make scree
  sc_objs <- original_and_robust(x,y)
  scr <- sc_objs$scree_ori
  if(out.rm) scr <- sc_objs$scree_rob

  #make alpha hull
  #d_ahull <- alphahull::ahull(x, y, a=alpha)
  ahull <- alphahull::ahull(scr$del, alpha=scr$alpha) # ahull_alpha) leave it out for now
  #d_ahull <- tibble::as_tibble(ahull$ashape.obj$edges)

  d_ahull <- tibble::tibble(x1 = ahull[["xahull"]][,1][ahull[["arcs"]][,7]],
                           y1 = ahull[["xahull"]][,2][ahull[["arcs"]][,7]],
                           x2 = ahull[["xahull"]][,1][ahull[["arcs"]][,8]],
                           y2 = ahull[["xahull"]][,2][ahull[["arcs"]][,8]])
  #scaled tibble
  d <- tibble::tibble(x=ahull[["xahull"]][,1], y=ahull[["xahull"]][,2])
  #d <- tibble::tibble(x=scr$del$x[,1], y=scr$del$x[,2])
  # d <- tibble::tibble(x=x, y=y)

  #make plot
  if(!fill){
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(data=d, ggplot2::aes(x, y),
                 colour = "black", alpha=alpha)+
      ggplot2::geom_segment(data=d_ahull,
                   ggplot2::aes(x=x1, xend=x2, y=y1, yend=y2))
  }
  if (fill){ #draws points on top of fill
    p <- ggplot2::ggplot() +
      ggplot2::geom_polygon(data=d_ahull,
                                   ggplot2::aes(x=x1, y=y1), alpha = 0.5) +
      ggplot2::geom_point(data=d, ggplot2::aes(x, y),
                          colour = "black", alpha=alpha)+
      ggplot2::geom_segment(data=d_ahull,
                            ggplot2::aes(x=x1, xend=x2, y=y1, yend=y2))
  }
  p
}

#' Drawing the MST
#'
#' This function will draw the MST for a
#' scatterplot.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param binner an optional function that bins the x and y vectors prior
#' to triangulation
#'  Can be:
#'   \itemize{
#'     \item `NULL`: no binning (use raw points)
#'     \item `"hex"` : hexagonal binning following the procedure in the
#'     graph-theoretic scagnostics paper (start 40x40, halve
#'              until <= 250 nonempty cells)
#'     \item a function: user-defined binner
#'     }
#' @param alpha character, numeric, or function. Controls the alpha radius.
#'   Valid character values are:
#'   \itemize{
#'     \item "rahman" (default): Rahman's MST-based middle-50% alpha
#'     \item "q90": 90th percentile of MST edge lengths
#'     \item "omega": graph-theoretic scagnostics alpha
#'   Alternatively:
#'     \item a numeric value giving a fixed alpha
#'     \item a function with no arguments that returns a single numeric alpha
#'     }
#' @param outlier_rm logical; if TRUE, iteratively trim large MST edges
#' @return A "gg" object that draws the plot's MST.
#' @examples
#' require(dplyr)
#' require(ggplot2)
#' data("features")
#' nl <- features %>% filter(feature == "positive")
#' draw_mst(nl$x, nl$y)
#' @export
draw_mst <- function(x, y, outlier_rm = FALSE, binner = NULL){
  UseMethod("draw_mst")
}

#' @rdname draw_mst
#' @export
draw_mst.default <- function(x, y, outlier_rm = FALSE,
                             binner = NULL){
  # ind1 <- ind2 <- connected <- x1 <- x2 <- y1 <- y2 <- NULL
  sc <- scree(x, y, outlier_rm, binner)
  draw_mst.scree(sc)
}

#' @rdname draw_mst
#' @export
draw_mst.scree <- function(x, y=NULL, outlier_rm = FALSE,
                           binner = NULL){
  mst <- gen_mst(x$del, x$weights)
  draw_mst.igraph(mst, x)
}

#' @rdname draw_mst
#' @export
draw_mst.igraph <- function(x, y, outlier_rm = FALSE,
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
draw_mst_old <- function(x, y, out.rm=TRUE) {
  ind1 <- ind2 <- connected <- x1 <- x2 <- y1 <- y2 <- NULL
  #build scree
  scree_obj <- original_and_robust(x, y)
  scree <- scree_obj$scree_ori
  if(out.rm) scree <- scree_obj$scree_rob

  #get mst
  MST <- gen_mst(scree$del, scree$weights)
  MST <- scree_obj$mst_ori
  if(out.rm) MST <- scree_obj$mst_rob

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
#' @param fill Fill the polygon
#' @param binner an optional function that bins the x and y vectors prior
#' to triangulation
#'  Can be:
#'   \itemize{
#'     \item `NULL`: no binning (use raw points)
#'     \item `"hex"` : hexagonal binning following the procedure in the
#'     graph-theoretic scagnostics paper (start 40x40, halve
#'              until <= 250 nonempty cells)
#'     \item a function: user-defined binner
#'     }
#' @param alpha character, numeric, or function. Controls the alpha radius.
#'   Valid character values are:
#'   \itemize{
#'     \item "rahman" (default): Rahman's MST-based middle-50% alpha
#'     \item "q90": 90th percentile of MST edge lengths
#'     \item "omega": graph-theoretic scagnostics alpha
#'   Alternatively:
#'     \item a numeric value giving a fixed alpha
#'     \item a function with no arguments that returns a single numeric alpha
#'     }
#' @param outlier_rm logical; if TRUE, iteratively trim large MST edges
#' @return A "gg" object that draws the plot's convex hull.
#' @examples
#' require(dplyr)
#' require(ggplot2)
#' data("features")
#' nl <- features %>% filter(feature == "clusters")
#' draw_convexhull(nl$x, nl$y, fill=TRUE, outlier_rm = FALSE)
#' @export
draw_convexhull <- function(x, y, outlier_rm = FALSE, binner = NULL,
                            fill = FALSE){
  UseMethod("draw_convexhull")
}

#' @rdname draw_convexhull
#' @export
draw_convexhull.default <- function(x, y, outlier_rm = FALSE,
                                    binner = NULL, fill = FALSE){
  # x1 <- x2 <- y1 <- y2 <- NULL
  sc <- scree(x, y, outlier_rm, binner)
  draw_convexhull.scree(sc, fill = fill)
}

#' @rdname draw_convexhull
#' @export
draw_convexhull.scree <- function(x, y=NULL, outlier_rm = FALSE,
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

  p
  }
  }

#' @export
draw_convexhull_old <- function(x, y, fill = FALSE, out.rm=TRUE) {
  x1 <- x2 <- y1 <- y2 <- NULL # happy cran checks
  # make scree and convex hull
  scree_obj <- original_and_robust(x, y)
  sc <- scree_obj$scree_ori
  if(out.rm) sc <- scree_obj$scree_rob
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

  if (fill) # re-draws points on top of fill
    p <- p + ggplot2::geom_polygon(data=d_ends, alpha=0.5,
                                   ggplot2::aes(x=x1, y=y1)) +
         ggplot2::geom_segment(data = d_ends,
                          ggplot2::aes(x=x1, y=y1,
                                       xend=x2, yend=y2)) +
         ggplot2::geom_point(data = d,
                             ggplot2::aes(x, y))

  p
}
