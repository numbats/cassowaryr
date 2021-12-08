#' Compute angle adjusted striated measure using MST
#'
#' @param x numeric vector of x values, or an MST object
#' @param y numeric vector of y values, or a scree object
#' @return A "numeric" object that gives the plot's striated2 score.
#'
#' @examples
#'   require(ggplot2)
#'   require(dplyr)
#'   ggplot(features, aes(x=x, y=y)) +
#'      geom_point() +
#'      facet_wrap(~feature, ncol = 5, scales = "free")
#'   features %>% group_by(feature) %>% summarise(striated = sc_striated2(x,y))
#'   sc_striated2(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'
#' @export
sc_striated2 <- function(x, y) UseMethod("sc_striated2")

#' @rdname sc_striated2
#' @export
sc_striated2.default <- function(x, y){
  sc <- scree(x, y)
  sc_striated2.scree(sc)
}

#' @rdname sc_striated2
#' @export
sc_striated2.scree <- function(x, y = NULL) {
  stopifnot(is.null(y))
  y <- gen_mst(x$del, x$weights)
  sc_striated2.igraph(y, x)

}

#' @rdname sc_striated2
#' @export
sc_striated2.igraph <- function(x, y){
  vertex_counts <- igraph::degree(x)
  angs <- which(vertex_counts>=2)
  stri=0
  for(i in seq_len(length(angs))){
    adjs <- which(x[angs[i]]>0)
    points <- y$del$x[adjs,]
    origin <- y$del$x[angs[i],]
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

#' Compute adjusted clumpy measure using MST
#'
#' @param x numeric vector of x values
#' @param y numeric vector of y values
#' @return A "numeric" object that gives the plot's clumpy2 score.
#'
#' @examples
#'   require(ggplot2)
#'   require(dplyr)
#'   ggplot(features, aes(x=x, y=y)) +
#'      geom_point() +
#'      facet_wrap(~feature, ncol = 5, scales = "free")
#'   features %>% group_by(feature) %>% summarise(clumpy = sc_clumpy2(x,y))
#'   sc_clumpy2(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
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
sc_clumpy2.scree <- function(x, y=NULL) {
  mst <- gen_mst(x$del, x$weights)
  sc_clumpy2.igraph(mst,x)
}

#' @rdname sc_clumpy2
#' @export
sc_clumpy2.igraph <- function(x, y){
  #set stringy penalty
  vertex_counts <- igraph::degree(x)
  #technically stringy calc
  stringy <- sum(vertex_counts == 2) / (length(vertex_counts) - sum(vertex_counts == 1))
  stringy_pen <- ifelse(stringy>0.95, (1-stringy), 1)
  # get lower triangular matrix
  mstmat <- twomstmat(x, y)$lowertri

  #get the index of all the edges in the mst
  matind <- which(mstmat>0) # in whole matrix
  rows <- matind %% length(mstmat[,1]) # respective row
  cols <- (matind-1) %/% length(mstmat[,1]) +1 #respective col
  #for later calculations
  n_total <- length(matind) #get total number of edges for weighting later

  # identify big edges
  edge_order <- order(mstmat[matind]) #get edge order
  ind <- which.min(diff(sort(mstmat[matind], decreasing = TRUE))) #index of maximum difference
  big_ei <- which(mstmat[matind] %in% utils::head(sort(mstmat, decreasing=TRUE), ind)) #index in original mst of big edges

  #only keep index, rows and cols of between cluster indexs
  matind <- matind[big_ei]
  rows <- rows[big_ei]
  cols <- cols[big_ei]

  #place holder for clumpy values
  clumpy <- rep(0,length(matind)) #1 for each between cluster edge

  #get big edge weights for later calculation
  big_ew <- mstmat[matind]
  #set mst we are going to remove all the edges from
  mst_ej <- mstmat
  mst_ej[matind] = 0 #remove big edges

  # finc clumpy value for each big edge
  for(j in seq(length(rows))){

    #have two clusters sprawling out from the deleted edge (inex i)
    c1rowcol <- rows[j]
    c2rowcol <- cols[j]

    #remake index objects to edit within iteration
    matind_ej <- which(mst_ej>0)
    rows_ej <- matind_ej %% length(mst_ej[,1])
    cols_ej <- (matind_ej-1) %/% length(mst_ej[,1]) +1

    #initialise variable that checks if clusters have changed
    whilecheck <- c(c1rowcol,c2rowcol)
    quitloop = 0

    while(quitloop==0){

      #find matches in rows/columns to join connected vertices
      c1rowcol <- unique(c(c1rowcol,
                           rows_ej[which(cols_ej%in%c1rowcol)],
                           cols_ej[which(rows_ej%in%c1rowcol)]))
      c2rowcol <- unique(c(c2rowcol,
                           rows_ej[which(cols_ej%in%c2rowcol)],
                           cols_ej[which(rows_ej%in%c2rowcol)]))

      #check if indices are done updating
      if(setequal(whilecheck, c(c1rowcol,c2rowcol))){
        quitloop=1
      }

      #update while loop check
      whilecheck = c(c1rowcol,c2rowcol)

    }
    #get indicies in the matrix
    c1ind <- unique(c(matind_ej[which(cols_ej%in%c1rowcol)],
                      matind_ej[which(rows_ej%in%c1rowcol)]))
    c2ind <- unique(c(matind_ej[which(cols_ej%in%c2rowcol)],
                      matind_ej[which(rows_ej%in%c2rowcol)]))

    #pull weights from the matrix
    c1weights <- mst_ej[c1ind]
    c2weights <- mst_ej[c2ind]

    # get the number of points in each cluster
    len_c1 <- length(c1weights)
    len_c2 <- length(c2weights)
    # take the value of the smaller cluster (points wise)
    c_length <- ifelse(len_c1 > len_c2, len_c2, len_c1)
    # and get short edge
    short_edge <- ifelse(len_c1 > len_c2, stats::median(c2weights), stats::median(c1weights))
    # possible code for dealing with irregular clusters
    # short_edge <- ifelse(len_c1 == len_c2, max(c(c2weights,c1weights)), short_edge)

    # calculate clumpy value w penalty for uneven clusters
    uneven_pen <- sqrt((2*c_length)/(len_c1+len_c2))
    clumpy[j] <- stringy_pen*uneven_pen*(big_ew[j]/short_edge)
    #just setting this now will fix to something more appropriate later
    clumpy[j] <- ifelse(is.na(clumpy[j]), 1, clumpy[j]) #return 1 if all clusters are of size 1
  }
  #threshold to be considered clumpy is above 1
  value <- ifelse(mean(clumpy)< 1, 1, mean(clumpy)) #value <- mean(clumpy)+1
  #stringy penalty only for high values
  #return clumpy
  1-(1/value)
}


#' Compute  adjusted sparse measure using the alpha hull
#'
#' @param x numeric vector of x values
#' @param y numeric vector of y values
#' @return A "numeric" object that gives the plot's sparse2 score.
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
sc_sparse2.scree <- function(x, y=NULL) {
  ahull <- gen_alpha_hull(x$del, x$alpha)
  sc_sparse2.list(ahull)
}

#' @rdname sc_sparse2
#' @export
sc_sparse2.list <- function(x, y=NULL){
  if (x$length > 0)
    ahull_area <- alphahull::areaahull(x)
  else
    ahull_area <- 0
  1- ahull_area
}
