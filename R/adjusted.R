#' Compute the grid scanostic measure using MST
#'
#' The grid scagnsotic as defined in Adam Rahman's PhD thesis (2018).
#' The scagnostic identifies grid-like structures by counting the number of
#' 90 and 180 degree angles in the MST. This measure can be used as an
#' effective alternative to striated when computing scagnostics without binning.
#'
#' @inheritParams scree
#' @param epsilon the error tolerance allowed when deciding if the MST angles
#' are at a right angle or not
#' @return A numeric object that gives the plot's grid score.
#' @examples
#'   require(ggplot2)
#'   require(dplyr)
#'   ggplot(features, aes(x=x, y=y)) +
#'      geom_point() +
#'      facet_wrap(~feature, ncol = 5, scales = "free")
#'   features |> group_by(feature) |>
#'     summarise(grid1 = sc_grid(x,y),
#'               grid2 = sc_grid(x,y, epsilon=0.05))
#'   sc_grid(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'
#' @export
sc_grid <- function(x, y, epsilon,
                    out.rm = TRUE, binner =  "hex") UseMethod("sc_grid")


#' @export
sc_grid.default <- function(x, y, epsilon=0.01, out.rm = TRUE, binner =  "hex"){
  sc <- scree(x, y, out.rm = out.rm, binner = binner)
  sc_grid.scree(sc, y=NULL, epsilon)
}


#' @export
sc_grid.scree <- function(x, y, epsilon=0.01, out.rm = TRUE, binner =  "hex") {
  stopifnot(is.null(y))
  y <- gen_mst(x$del, x$weights)
  sc_grid.igraph(y, x, epsilon)

}


#' @export
sc_grid.igraph <- function(x, y, epsilon=0.01, out.rm = TRUE, binner =  "hex"){
  vertex_counts <- igraph::degree(x)
  angs <- which(vertex_counts>=2)
  grd=0
  for(i in seq_len(length(angs))){
    adjs <- which(x[angs[i]]>0)
    points <- y$del$x[adjs,]
    origin <- y$del$x[angs[i],]
    vects <- t(t(points)-origin)
    b =0
    for(j in seq(length(vects[,1])-1)){
      costheta <- (vects[j,]%*%vects[j+1,])/(sqrt(sum(vects[j,]^2))*sqrt(sum(vects[j+1,]^2)))
      b <- ifelse(any(c(costheta<(-1+epsilon), abs(costheta)<epsilon)), b+1, b)
    }
    grd <- grd + b
  }
  grd/(0.5*sum(vertex_counts) - 1)
}

#' Compute adjusted clumpy measure using MST
#'
#' This measure is defined in the cassowaryr paper by Mason, et al. (2025).
#' It is an alternative measure for clumpiness. It is the ratio of the
#' between cluster edges and the within cluster edges. It is a good
#' alternative measure to clumpy when binning is removed as a pre-processing
#' step.
#'
#' @inheritParams scree
#' @return A numeric object that gives the plot's adjusted clumpy score.
#' @examples
#'   require(ggplot2)
#'   require(dplyr)
#'   ggplot(features, aes(x=x, y=y)) +
#'      geom_point() +
#'      facet_wrap(~feature, ncol = 5, scales = "free")
#'   features |> group_by(feature) |> summarise(clumpy = sc_clumpy2(x,y))
#'   sc_clumpy2(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'
#' data <- features |> filter(feature == "clusters")
#' x <- data$x
#' y <- data$y
#'
#' # plot it
#' ggplot() +
#'   geom_point(aes(x = x, y = y))
#'
#' # calculate using vectors
#' sc_clumpy2(x, y)
#'
#' @export
sc_clumpy2 <- function(x, y, out.rm = TRUE, binner =  "hex") UseMethod("sc_clumpy2")


#' @export
sc_clumpy2.default <- function(x, y, out.rm = TRUE, binner =  "hex"){
  sc <- scree(x, y, out.rm = out.rm, binner = binner)
  sc_clumpy2.scree(sc)
}


#' @export
sc_clumpy2.scree <- function(x, y=NULL, out.rm = FALSE, binner = NULL) {
  mst <- gen_mst(x$del, x$weights)
  sc_clumpy2.igraph(mst,x)
}


#' @export
sc_clumpy2.igraph <- function(x, y=NULL, out.rm = TRUE, binner =  "hex"){
  # Get...
  # edge weights
  mst_weights <- igraph::E(x)$weight
  # number of edges
  n_total <- length(mst_weights)
  # Lower triangular matrix
  mstmat <- lower.tri(matrix(x[], nrow = n_total + 1))

  #get the index of all the edges in the mst
  matind <- which(mstmat>0) # in whole matrix
  rows <- matind %% length(mstmat[,1]) # respective row
  cols <- (matind-1) %/% length(mstmat[,1]) +1 #respective col
  print(matind)
  print(rows)
  print(cols)

  # Sort and get length difference
  #index of maximum difference
  ind <- which.min(diff(sort(mst_weights, decreasing = TRUE)))

  print(ind)
  #index in original mst of big edges
  big_ei <- which(mst_weights %in% utils::head(sort(mst_weights, decreasing=TRUE), ind))
  print(big_ei)
  #only keep index, rows and cols of between cluster indexs
  matind <- matind[big_ei]
  rows <- rows[big_ei]
  cols <- cols[big_ei]
  print(matind)
  print(rows)
  print(cols)
  #place holder for clumpy values
  clumpy <- rep(0,length(matind)) #1 for each between cluster edge
  print(clumpy)
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
    clumpy[j] <- uneven_pen*(big_ew[j]/short_edge)
    #return 1 if all clusters are of size 1
    clumpy[j] <- ifelse(is.na(clumpy[j]), 1, clumpy[j])
  }

  #threshold to be considered clumpy is above 1
  value <- ifelse(mean(clumpy)< 1, 1, mean(clumpy))

  #return clumpy
  1-(1/value)
}


#' Compute  adjusted sparse measure using the alpha hull
#'
#' The sparse2 measure created for cassowaryr
#' The measure calculates the sparsity of the plot as 1-area(ahull).
#'
#' @inheritParams scree
#' @return A numeric object that gives the plot's adjusted sparse score.
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
sc_sparse2 <- function(x, y, out.rm = TRUE, binner =  "hex") UseMethod("sc_sparse2")


#' @export
sc_sparse2.default <- function(x, y, out.rm = TRUE, binner =  "hex"){
  sc <- scree(x, y, out.rm = out.rm, binner = binner)
  sc_sparse2.scree(sc)
}


#' @export
sc_sparse2.scree <- function(x, y=NULL, out.rm = TRUE, binner =  "hex") {
  ahull <- gen_alpha_hull(x$del, x$alpha)
  sc_sparse2.list(ahull)
}


#' @export
sc_sparse2.list <- function(x, y=NULL, out.rm = TRUE, binner =  "hex"){
  if (x$length > 0)
    ahull_area <- alphahull::areaahull(x)
  else
    ahull_area <- 0
  1- ahull_area
}


