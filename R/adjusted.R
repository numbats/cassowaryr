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
  # get lower triangular matrix
  mstmat <- twomstmat(mst, sc)$lowertri

  #get the index of all the edges in the mst
  matind <- which(mstmat>0) # in whole matrix
  rows <- matind %% length(mstmat[,1]) # respective row
  cols <- (matind-1) %/% length(mstmat[,1]) +1 #respective col
  #for later calculations
  n_total <- length(matind) #get total number of edges for weighting later

  # identify big edges
  edge_order <- order(mstmat[matind]) #get edge order
  ind <- which.min(diff(sort(mstmat[matind], decreasing = TRUE))) #index of maximum difference
  big_ei <- which(mstmat[matind] %in% head(sort(mstmat, decreasing=TRUE), ind)) #index in original mst of big edges

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
    short_edge <- ifelse(len_c1 > len_c2, median(c2weights), median(c1weights))
    # possible code for dealing with irregular clusters
    # short_edge <- ifelse(len_c1 == len_c2, max(c(c2weights,c1weights)), short_edge)

    # calculate clumpy value w penalty for uneven clusters
    uneven <- ((2*c_length)/(len_c1+len_c2))
    clumpy[j] <- uneven*(big_ew[j]/short_edge)
  }
  mean(clumpy)
}

sc_clumpy.igraph <- function(mymst, x){

  #lower triangular matrix
  mstmat <- twomstmat(mymst,x)$lowertri

  #make index variables to iterate through
  matind <- which(mstmat>0)
  rows <- matind %% length(mstmat[,1])
  cols <- (matind-1) %/% length(mstmat[,1]) +1
  clumpy <- rep(0,length(matind))

  for(j in seq(length(rows))){
    #set mst we are going to remove all the edges from
    mst_ej <- mstmat

    #have two clusters sprawling out from the deleted edge (inex i)
    c1rowcol <- rows[j]
    c2rowcol <- cols[j]

    #get weight of ej
    ej_weight <- mst_ej[matind[j]]

    #remove ej and all values in mst that are greater than ej
    mst_ej[which(mst_ej>=ej_weight)] = 0

    #remake index objects to edit within iteration
    matind_ej <- which(mst_ej>0)
    rows_ej <- matind_ej %% length(mst_ej[,1])
    cols_ej <- (matind_ej-1) %/% length(mst_ej[,1]) +1

    #initialise variable that checks if clusters have changed
    whilecheck <- c(c1rowcol,c2rowcol)
    quitloop = 0

    #add in new indices until both clusters are full
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

    #get matrix indices of each of the clusters
    c1ind <- unique(c(matind_ej[which(cols_ej%in%c1rowcol)],
                      matind_ej[which(rows_ej%in%c1rowcol)]))
    c2ind <- unique(c(matind_ej[which(cols_ej%in%c2rowcol)],
                      matind_ej[which(rows_ej%in%c2rowcol)]))

    #get edge weights for each cluster
    c1weights <- mst_ej[c1ind]
    c2weights <- mst_ej[c2ind]

    #set K weight value
    ek <- max(c(0,c1weights)) #max(c1weights, na.rm=TRUE)
    em <- max(c(0,c2weights)) #max(c2weights, na.rm=TRUE)
    ek_weight <- ifelse(length(c1weights)<length(c2weights), ek, em)


    #calculate this clumpy value
    clumpy[j] <- ifelse(ek_weight==0, 0, 1 - (ek_weight/ej_weight)) # ifelse(is.finite(1 - ek_weight/ej_weight), 1 - ek_weight/ej_weight, 0)
  }

  #remove NA and return final clumpy measure
  max(clumpy, na.rm=TRUE)

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
