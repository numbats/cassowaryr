
#' Compute stringy scagnostic measure using MST
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_stringy(anscombe$x1, anscombe$y1)
#'   sc_stringy(anscombe$x2, anscombe$y2)
#'   sc_stringy(anscombe$x3, anscombe$y3)
#'   sc_stringy(anscombe$x4, anscombe$y4)
#'
#' @export
sc_stringy <- function(x, y) UseMethod("sc_stringy")

#' @rdname sc_stringy
#' @export
sc_stringy.scree <- function(x, y = NULL) {
  stopifnot(is.null(y))
  mst <- gen_mst(x$del, x$weights)
  sc_stringy.igraph(mst)
}

#' @rdname sc_stringy
#' @export
sc_stringy.default <- function(x, y){
  sc <- scree(x, y)
  sc_stringy.scree(sc)
}

#' @rdname sc_stringy
#' @export
sc_stringy.igraph <- function(mst){
  vertex_counts <- igraph::degree(mst)
  sum(vertex_counts == 2) / (length(vertex_counts) - sum(vertex_counts == 1))
}

#' Compute stirated scagnostic measure using MST
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_striated(anscombe$x1, anscombe$y1)
#'   sc_striated(anscombe$x2, anscombe$y2)
#'   sc_striated(anscombe$x3, anscombe$y3)
#'   sc_striated(anscombe$x4, anscombe$y4)
#'
#'   ggplot(datasaurus_dozen, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~dataset, ncol=3, scales = "free")
#'   sc_striated(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'   sc_striated(datasaurus_dozen_wide$dino_x, datasaurus_dozen_wide$dino_y)
#'   sc_striated(datasaurus_dozen_wide$dots_x, datasaurus_dozen_wide$dots_y)
#'   sc_striated(datasaurus_dozen_wide$h_lines_x, datasaurus_dozen_wide$h_lines_y)
#'
#' @export
sc_striated <- function(x, y) UseMethod("sc_striated")

#' @rdname sc_striated
#' @export
sc_striated.scree <- function(x, y = NULL) {
  mst <- gen_mst(x$del, x$weights)
  sc_striated.igraph(mst, x)

  }

#' @rdname sc_striated
#' @export
sc_striated.default <- function(x, y){
  sc <- scree(x, y)
  sc_striated.scree(sc)
}

#' @rdname sc_striated
#' @export
sc_striated.igraph <- function(mst, x){
  vertex_counts <- igraph::degree(mst)
  angs <- which(vertex_counts==2)
  angles_vect <- numeric(length(angs))
  for(i in seq_len(length(angs))){
    adjs <- which(mst[angs[i]]>0)
    points <- x$del$x[adjs,]
    origin <- x$del$x[angs[i],]
    vects <- t(t(points)-origin)
    angles_vect[i] <- (vects[1,]%*%vects[2,])/(prod(mst[angs[i]][adjs]))
  }
  (sum(ifelse(angles_vect<(-0.75),1,0)))/length(vertex_counts)
}

#' Compute angle adjusted stirated measure using MST
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_striated_adjusted(anscombe$x1, anscombe$y1)
#' @export
sc_striated_adjusted <- function(x, y) UseMethod("sc_striated_adjusted")

#' @rdname sc_striated_adjusted
#' @export
sc_striated_adjusted.scree <- function(x, y = NULL) {
  mst <- gen_mst(x$del, x$weights)
  sc_striated_adjusted.igraph(mst, x)

}

#' @rdname sc_striated_adjusted
#' @export
sc_striated_adjusted.default <- function(x, y){
  sc <- scree(x, y)
  sc_striated_adjusted.scree(sc)
}

#' @rdname sc_striated_adjusted
#' @export
sc_striated_adjusted.igraph <- function(mst, x){
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

#' Compute clumpy scagnostic measure using MST
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_clumpy(anscombe$x1, anscombe$y1)
#'   sc_clumpy(anscombe$x2, anscombe$y2)
#'   sc_clumpy(anscombe$x3, anscombe$y3)
#'   sc_clumpy(anscombe$x4, anscombe$y4)
#'
#'   ggplot(datasaurus_dozen, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~dataset, ncol=3, scales = "free")
#'   sc_clumpy(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'   sc_clumpy(datasaurus_dozen_wide$dino_x, datasaurus_dozen_wide$dino_y)
#'   sc_clumpy(datasaurus_dozen_wide$dots_x, datasaurus_dozen_wide$dots_y)
#'   sc_clumpy(datasaurus_dozen_wide$h_lines_x, datasaurus_dozen_wide$h_lines_y)
#'
#' @export
sc_clumpy <- function(x, y) UseMethod("sc_clumpy")

#' @rdname sc_clumpy
#' @export
sc_clumpy.scree <- function(x, y = NULL) {
  mymst <- gen_mst(x$del, x$weights)
  sc_clumpy.igraph(mymst,x)
}

#' @rdname sc_clumpy
#' @export
sc_clumpy.default <- function(x, y){
  sc <- scree(x, y)
  sc_clumpy.scree(sc)
}

#' @rdname sc_clumpy
#' @export
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

#' @export
sc_clumpy_adjusted <- function(x, y) UseMethod("sc_clumpy_adjusted")

#' @rdname sc_clumpy_adjusted
#' @export
sc_clumpy_adjusted.default <- function(x, y){
  sc <- scree(x, y)
  sc_clumpy_adjusted.scree(sc)
}

#' @rdname sc_clumpy_adjusted
#' @export
sc_clumpy_adjusted.scree <- function(x, y = NULL) {
  mymst <- gen_mst(x$del, x$weights)
  sc_clumpy_adjusted.igraph(mymst,x)
}

#' @rdname sc_clumpy_adjusted
#' @export
sc_clumpy_adjusted.igraph <- function(mst, sc){
  mstmat <- twomstmat(mst, sc)$lowertri
  edges <- sort(mstmat[which(mstmat>0)])
  ind <- which.max(diff(edges))
  bigedges <- mean(edges[(ind+1):length(edges)])
  smalledges <- max(edges[1:ind])
  1- smalledges/bigedges
}


#' Compute sparse scagnostic measure using MST
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_sparse(anscombe$x1, anscombe$y1)
#'   sc_sparse(anscombe$x2, anscombe$y2)
#'   sc_sparse(anscombe$x3, anscombe$y3)
#'   sc_sparse(anscombe$x4, anscombe$y4)
#'   ggplot(datasaurus_dozen, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~dataset, ncol=3, scales = "free")
#'   sc_sparse(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'   sc_sparse(datasaurus_dozen_wide$bullseye_x, datasaurus_dozen_wide$bullseye_y)
#'   sc_sparse(datasaurus_dozen_wide$circle_x, datasaurus_dozen_wide$circle_y)
#'   sc_sparse(datasaurus_dozen_wide$dino_x, datasaurus_dozen_wide$dino_y)
#'   sc_sparse(datasaurus_dozen_wide$dots_x, datasaurus_dozen_wide$dots_y)
#'   sc_sparse(datasaurus_dozen_wide$h_lines_x, datasaurus_dozen_wide$h_lines_y)
#'   sc_sparse(datasaurus_dozen_wide$high_lines_x, datasaurus_dozen_wide$high_lines_y)
#'   sc_sparse(datasaurus_dozen_wide$slant_down_x, datasaurus_dozen_wide$slant_up_y)
#'   sc_sparse(datasaurus_dozen_wide$star_x, datasaurus_dozen_wide$star_y)
#'   sc_sparse(datasaurus_dozen_wide$v_lines_x, datasaurus_dozen_wide$v_lines_y)
#'   sc_sparse(datasaurus_dozen_wide$wide_lines_x, datasaurus_dozen_wide$wide_lines_y)
#'   sc_sparse(datasaurus_dozen_wide$x_shape_x, datasaurus_dozen_wide$x_shape_y)
#'
#' @export
sc_sparse <- function(x, y) UseMethod("sc_sparse")

#' @rdname sc_sparse
#' @export
sc_sparse.scree <- function(x, y = NULL) {
  #generate vector of MST edges
  mymst <- gen_mst(x$del, x$weights)
  sc_sparse.igraph(mymst,x)
}

#' @rdname sc_sparse
#' @export
sc_sparse.default <- function(x, y){
  sc <- scree(x, y)
  sc_sparse.scree(sc)
}

#' @rdname sc_sparse
#' @export
sc_sparse.igraph <- function(mymst, x){
  mstmat <- twomstmat(mymst,x)$lowertri
  edges <- mstmat[which(mstmat>0)]
  #calculate sparse value
  sort(edges)[floor(0.9*length( edges))]
}


#' Compute skewed scagnostic measure using MST
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_skewed(anscombe$x1, anscombe$y1)
#'   sc_skewed(anscombe$x2, anscombe$y2)
#'   sc_skewed(anscombe$x3, anscombe$y3)
#'   sc_skewed(anscombe$x4, anscombe$y4)
#'   ggplot(datasaurus_dozen, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~dataset, ncol=3, scales = "free")
#'   sc_skewed(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'   sc_skewed(datasaurus_dozen_wide$bullseye_x, datasaurus_dozen_wide$bullseye_y)
#'   sc_skewed(datasaurus_dozen_wide$circle_x, datasaurus_dozen_wide$circle_y)
#'   sc_skewed(datasaurus_dozen_wide$dino_x, datasaurus_dozen_wide$dino_y)
#'   sc_skewed(datasaurus_dozen_wide$dots_x, datasaurus_dozen_wide$dots_y)
#'   sc_skewed(datasaurus_dozen_wide$h_lines_x, datasaurus_dozen_wide$h_lines_y)
#'   sc_skewed(datasaurus_dozen_wide$high_lines_x, datasaurus_dozen_wide$high_lines_y)
#'   sc_skewed(datasaurus_dozen_wide$slant_down_x, datasaurus_dozen_wide$slant_up_y)
#'   sc_skewed(datasaurus_dozen_wide$star_x, datasaurus_dozen_wide$star_y)
#'   sc_skewed(datasaurus_dozen_wide$v_lines_x, datasaurus_dozen_wide$v_lines_y)
#'   sc_skewed(datasaurus_dozen_wide$wide_lines_x, datasaurus_dozen_wide$wide_lines_y)
#'   sc_skewed(datasaurus_dozen_wide$x_shape_x, datasaurus_dozen_wide$x_shape_y)
#'
#' @export
sc_skewed <- function(x, y) UseMethod("sc_skewed")

#' @rdname sc_skewed
#' @export
sc_skewed.scree <- function(x, y = NULL) {
  #generate vector of MST edges
  mymst <- gen_mst(x$del, x$weights)
  sc_skewed.igraph(mymst, x)
}

#' @rdname sc_skewed
#' @export
sc_skewed.default <- function(x, y){
  sc <- scree(x, y)
  sc_skewed.scree(sc)
}

#' @rdname sc_skewed
#' @export
sc_skewed.igraph <- function(mymst, x){
  mstmat <- twomstmat(mymst,x)$lowertri
  edges <- mstmat[which(mstmat>0)]

  # find quantiles
  q10 <- sort(edges)[floor(0.1*length( edges))]
  q50 <- sort(edges)[floor(0.5*length( edges))]
  q90 <- sort(edges)[floor(0.9*length( edges))]

  # calculate skewed value
  (q90-q50)/(q90-q10)

}


#' Compute outlying scagnostic measure using MST
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_outlying(anscombe$x1, anscombe$y1)
#'   sc_outlying(anscombe$x2, anscombe$y2)
#'   sc_outlying(anscombe$x3, anscombe$y3)
#'   sc_outlying(anscombe$x4, anscombe$y4)
#'   ggplot(datasaurus_dozen, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~dataset, ncol=3, scales = "free")
#'   sc_outlying(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'   sc_outlying(datasaurus_dozen_wide$bullseye_x, datasaurus_dozen_wide$bullseye_y)
#'   sc_outlying(datasaurus_dozen_wide$circle_x, datasaurus_dozen_wide$circle_y)
#'   sc_outlying(datasaurus_dozen_wide$dino_x, datasaurus_dozen_wide$dino_y)
#'   sc_outlying(datasaurus_dozen_wide$dots_x, datasaurus_dozen_wide$dots_y)
#'   sc_outlying(datasaurus_dozen_wide$h_lines_x, datasaurus_dozen_wide$h_lines_y)
#'   sc_outlying(datasaurus_dozen_wide$high_lines_x, datasaurus_dozen_wide$high_lines_y)
#'   sc_outlying(datasaurus_dozen_wide$slant_down_x, datasaurus_dozen_wide$slant_up_y)
#'   sc_outlying(datasaurus_dozen_wide$star_x, datasaurus_dozen_wide$star_y)
#'   sc_outlying(datasaurus_dozen_wide$v_lines_x, datasaurus_dozen_wide$v_lines_y)
#'   sc_outlying(datasaurus_dozen_wide$wide_lines_x, datasaurus_dozen_wide$wide_lines_y)
#'   sc_outlying(datasaurus_dozen_wide$x_shape_x, datasaurus_dozen_wide$x_shape_y)
#'
#' @export
sc_outlying <- function(x, y) UseMethod("sc_outlying")

#' @rdname sc_outlying
#' @export
sc_outlying.scree <- function(x, y = NULL) {
  #generate vector of MST edges
  mymst <- gen_mst(x$del, x$weights)
  sc_outlying.igraph(mymst, x)
}

#' @rdname sc_outlying
#' @export
sc_outlying.default <- function(x, y){
  sc <- scree(x, y)
  sc_outlying.scree(sc)
}

#' @export
sc_outlying.igraph <- function(mymst, x){
  #input: original mst (mymst) and scree object (x)
  #output: outlying mst value

  #make into matrix
  twomst <- twomstmat(mymst,x)
  mstmat <- twomst$mat
  mstlowertri <- twomst$lowertri

  #identify outliers
  outliers <- outlying_identify(mymst, x)

  #calculate outlying value
  outlier_e <- sum(mstmat[outliers,]) #sum of edges of outlying points
  overlap <- sum(mstmat[as.matrix(expand.grid(outliers,outliers))]) #outerliers connected to outliers
  numer <- outlier_e - 0.5*overlap #overlap double counts
  denom <-  0.5*sum(mstmat)
  numer/denom
}

gen_mst <- function(del, weights) {
  edges <- del$mesh[, c("ind1", "ind2")]
  graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
  graph <- igraph::set_edge_attr(graph, "weight", value = weights)
  igraph::mst(graph, weights =  igraph::E(graph)$weight)
}

outlying_identify <- function(mst, scr){
  #input: takes a mst and scree
  #output: rown number of

  #get matrix and upper triangular matrix
  matlist <- twomstmat(mst,scr)
  mstmat <- matlist$mat
  mstmat_lt <- matlist$lowertri

  #calculate w value
  edges <- mstmat_lt[which(mstmat_lt>0)]
  w <- psi(edges)

  #set values above w to 0 in matrix to find outlying
  mstmat_check <- mstmat
  mstmat_check[mstmat>w]=0

  #row sum of matrix, if 0 all edges are above this value
  rowsum <- mstmat_check%*%rep(1, length(mstmat_check[1,]))
  #return the outlying observation
  which(rowsum==0)
}

original_and_robust <- function(x,y){
  #input: data for 2 variables x and y
  #output: list of scree and MST objects

  #construct original scree and MST
  sc_original <- scree(x,y)
  mst_original <- gen_mst(sc_original$del, sc_original$weights)

  #identify outliers
  outliers <- outlying_identify(mst_original, sc_original)

  #set outlier removed to original in case of no outliers
  sc_robust <- sc_original
  mst_robust <- mst_original

  #outlier removed scree and mst
  if(length(outliers)>0){
    new_x <- x[-outliers]
    new_y <- y[-outliers]
    if(sd(new_x) == 0 | sd(new_y) == 0) return(NULL)
    #recalculate scree and MST
    sc_robust <- scree(new_x,new_y)
    mst_robust <- gen_mst(sc_robust$del, sc_robust$weights)
  }

  #output 4 objects as a list
  structure(
    list(
      scree_ori = sc_original,
      mst_ori  = mst_original,
      scree_rob = sc_robust,
      mst_rob = mst_robust
    ))
}

twomstmat <- function(mst, scr){
  #input: mst and scree
  #output: mst full and lower triangular matrices

  #make into upper tri-matrix
  mst_mat <- matrix(mst[], nrow=length(scr[["del"]][["x"]][,1]))
  mst_uppertri <- mst_mat
  mst_uppertri[upper.tri(mst_mat, diag = FALSE)]=0

  #output matrix and upper triangular matrix
  structure(
    list(mat = mst_mat,
         lowertri = mst_uppertri))
}

