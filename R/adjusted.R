#' Compute angle adjusted stirated measure using MST
#' @param x numeric vector of x values
#' @param y numeric vector of y values
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_striated2(anscombe$x1, anscombe$y1)
#' @export
sc_striated2 <- function(x, y) UseMethod("sc_striated2")

#' @rdname sc_striated2
#' @export
sc_striated2.scree <- function(x, y = NULL) {
  mst <- gen_mst(x$del, x$weights)
  sc_striated2.igraph(mst, x)

}

#' @rdname sc_striated2
#' @export
sc_striated2.default <- function(x, y){
  sc <- scree(x, y)
  sc_striated2.scree(sc)
}

#' @rdname sc_striated2
#' @export
sc_striated2.igraph <- function(mst, x){
  vertex_counts <- igraph::degree(mst)
  angs <- which(vertex_counts>=2)
  stri=0
  for(i in seq_len(length(angs))){
    adjs <- which(mst[angs[i]]>0)
    points <- x$del$x[adjs,]
    origin <- x$del$x[angs[i],]
    vects <- t(t(points)-origin)
    b =0
    for(j in seq(length(vects[,1])-1)){
      costheta <- (vects[j,]%*%vects[j+1,])/(sqrt(sum(vects[j,]^2))*sqrt(sum(vects[j+1,]^2)))
      b <- ifelse(any(c(costheta<(-0.99), abs(costheta)<0.01)), b+1, b)
    }
    stri <- stri + b
  }
  stri/(0.5*sum(vertex_counts) - 1)
}

#' Compute angle adjusted clumpy measure using MST
#' @param x numeric vector of x values
#' @param y numeric vector of y values
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_clumpy2(anscombe$x1, anscombe$y1)
#'
#' @export
sc_clumpy2 <- function(x, y) UseMethod("sc_clumpy2")

#' @rdname sc_clumpy2
#' @export
sc_clumpy2.default <- function(x, y){
  sc <- scree(x, y)
  sc_clumpy2.scree(sc)
}

#' @rdname sc_clumpy2
#' @export
sc_clumpy2.scree <- function(x, y = NULL) {
  mymst <- gen_mst(x$del, x$weights)
  sc_clumpy2.igraph(mymst,x)
}

#' @rdname sc_clumpy2
#' @export
sc_clumpy2.igraph <- function(mst, sc){
  mstmat <- twomstmat(mst, sc)$lowertri
  edges <- sort(mstmat[which(mstmat>0)])
  ind <- which.max(diff(edges))
  bigedges <- mean(edges[(ind+1):length(edges)])
  smalledges <- max(edges[1:ind])
  1- smalledges/bigedges
}


#' Compute angle adjusted sparse measure using the Alpha Hull
#'
#' @param x numeric vector of x values
#' @param y numeric vector of y values
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_sparse2(anscombe$x1, anscombe$y1)
#'
#' @export
sc_sparse2 <- function(x, y) UseMethod("sc_sparse2")

#' @rdname sc_sparse2
#' @export
sc_sparse2.default <- function(x, y){
  sc <- scree(x, y)
  sc_sparse2.scree(sc)
}

#' @rdname sc_sparse2
#' @export
sc_sparse2.scree <- function(x, y = NULL) {
  stopifnot(is.null(y))
  ahull <- gen_alpha_hull(x$del, x$alpha)
  sc_sparse2.list(ahull)
}

#' @rdname sc_sparse2
#' @export
sc_sparse2.list <- function(ahull){
  if (ahull$length > 0)
    ahull_area <- alphahull::areaahull(ahull)
  else
    ahull_area <- 0
  1- ahull_area
}
