
#' Compute stringy scagnostic measure using MST
#'
#' This measure identifies a “stringy” shape with no branches, such as a thin
#' line of data. The stringy function is defined in Scagnostics Distributions by
#' Wilkinson & Wills (2008). It is calculated by comparing the number of
#' vertices of degree two in the MST with the total number of vertices
#' in the MST, dropping those of degree one. The stringy2 function
#' is defined in Graph Theoretic Scagnostics, Wilkinson (2005). It is
#' the legnth of the longest shortest path through the MST divided by the
#' sum of all edge lengths in the MST.

#'
#' @inheritParams scree
#' @return A numeric object that gives the plot's stringy score.
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
sc_stringy <- function(x, y, out.rm = TRUE, binner = "hex") UseMethod("sc_stringy")


#' @export
sc_stringy.default <- function(x, y, out.rm = TRUE, binner = "hex"){
  #input: x and y are vectors
  sc <- scree(x, y, out.rm = out.rm, binner = binner)
  if (is.null(sc$del)) {
    dc <- sc_dcor(x,y)
    if (dc >= 1 - 1e-8) {
      return(1)
    }

  }
  sc_stringy.scree(sc)
}


#' @export
sc_stringy.scree <- function(x, y, out.rm = TRUE, binner = "hex") {
  #input: x is a scree, no y
  mst <- gen_mst(x$del, x$weights)
  sc_stringy.igraph(mst)
}


#' @export
sc_stringy.igraph <- function(x, y, out.rm = TRUE, binner = "hex") {
  #input: x is the MST igraph object
  vertex_counts <- igraph::degree(x)
  sum(vertex_counts == 2) / (length(vertex_counts) - sum(vertex_counts == 1))
}

#' @rdname sc_stringy
#' @export
sc_stringy2 <- function(x, y, out.rm = TRUE, binner = "hex") UseMethod("sc_stringy2")


#' @export
sc_stringy2.default <- function(x, y, out.rm = TRUE, binner = "hex"){
  #input: x and y are vectors
  sc <- scree(x, y, out.rm = out.rm, binner = binner)
  if (is.null(sc$del)) {
    dc <- sc_dcor(x,y)
    if (dc >= 1 - 1e-8) {
      return(1)
    }

  }
  sc_stringy2.scree(sc)
}


#' @export
sc_stringy2.scree <- function(x, y=NULL, out.rm = FALSE, binner = NULL) {
  #input: x is a scree, no y
  mst <- gen_mst(x$del, x$weights)
  sc_stringy2.igraph(mst)
}


#' @export
sc_stringy2.igraph <- function(x, y=NULL, out.rm = FALSE, binner = NULL) {
  #input: x is the MST igraph object
  diameter <- igraph::get_diameter(x)
  length(diameter) / (length(x) - 1)
}

#' Compute striated scagnostic measure using MST
#'
#' This measure identifies features such as discreteness by finding parallel
#' lines. It was first defined in Graph Theoretic Scagnostics, Wilkinson,
#' et al. (2005). It is calculated by counting the proportion of vertices
#' with only two edges that have an inner angle approximately between
#' 135 and 220 degrees.
#'
#' @inheritParams scree
#' @return A numeric object that gives the plot's striated score.
#' @examples
#'   require(ggplot2)
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
#' @export
sc_striated <- function(x, y, out.rm = TRUE, binner = "hex") UseMethod("sc_striated")


#' @export
sc_striated.default <- function(x, y, out.rm = TRUE, binner = "hex"){
  sc <- scree(x, y, out.rm = out.rm, binner = binner)
  sc_striated.scree(sc)
}


#' @export
sc_striated.scree <- function(x, y=NULL, out.rm = FALSE, binner = NULL) {
  #input: x is scree
  mst <- gen_mst(x$del, x$weights)
  sc_striated.igraph(mst, x)
}



#' @export
sc_striated.igraph <- function(x, y, out.rm = TRUE, binner = "hex"){
  #input: x is MST, y is scree
  vertex_counts <- igraph::degree(x)
  angs <- which(vertex_counts==2)
  angles_vect <- numeric(length(angs))
  for(i in seq_len(length(angs))){
    adjs <- which(x[angs[i]]>0)
    points <- y$del$x[adjs,]
    origin <- y$del$x[angs[i],]
    vects <- t(t(points)-origin)
    angles_vect[i] <- (vects[1,]%*%vects[2,])/(prod(x[angs[i]][adjs]))
  }
  (sum(ifelse(angles_vect<(-0.75),1,0)))/length(vertex_counts)
}

#' Compute clumpy scagnostic measure using MST
#'
#' This measure is used to detect clustering and is calculated through an
#' iterative process. It was first defined in Graph Theoretic Scagnostics,
#' Wilkinson, et al. (2005). First an edge J is selected and removed from the
#' MST. From the two spanning trees that are created by this break, we select
#' the largest edge from the smaller tree (K). The length of this edge (K) is
#' compared to the removed edge (J) giving a clumpy measure for this edge.
#' This process is repeated for every edge in the MST and the final clumpy
#' measure is the maximum of this value over all edges.
#'
#' @inheritParams scree
#' @return A numeric object that gives the plot's clumpy score.
#' @examples
#'   require(ggplot2)
#'   require(dplyr)
#'   ggplot(features, aes(x=x, y=y)) +
#'      geom_point() +
#'      facet_wrap(~feature, ncol = 5, scales = "free")
#'   features |> group_by(feature) |> summarise(clumpy = sc_clumpy(x,y))
#'   sc_clumpy(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'
#' @export
sc_clumpy <- function(x, y, out.rm = TRUE, binner = "hex") UseMethod("sc_clumpy")


#' @export
sc_clumpy.default <- function(x, y, out.rm = TRUE, binner = "hex"){
  sc <- scree(x, y, out.rm = out.rm, binner = binner)
  sc_clumpy.scree(sc)
}


#' @export
sc_clumpy.scree <- function(x, y=NULL, out.rm = FALSE, binner = NULL) {
  #input: x is a scree
  mst <- gen_mst(x$del, x$weights)
  sc_clumpy.igraph(mst,x)
}



#' @export
sc_clumpy.igraph <- function(x, y, out.rm = TRUE, binner = "hex"){
  #input: x is the MST, y is the scree
  #lower triangular matrix
  mstmat <- twomstmat(x,y)$lowertri

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

#' Compute sparse scagnostic measure using MST
#'
#' Identifies if the data is confined to a small number of locations on the
#' plane. It was first defined in Scagnostics Distributions by
#' Wilkinson & Wills (2008). It is calculated as the 90th percentile
#' of MST edge lengths
#'
#' @inheritParams scree
#' @return A numeric object that gives the plot's sparse score.
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   ggplot(datasaurus_dozen, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~dataset, ncol=3, scales = "free")
#'   sc_sparse(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'   sc_sparse(datasaurus_dozen_wide$circle_x, datasaurus_dozen_wide$circle_y)
#'   sc_sparse(datasaurus_dozen_wide$dino_x, datasaurus_dozen_wide$dino_y)
#'
#' @export
sc_sparse <- function(x, y, out.rm = TRUE, binner = "hex") UseMethod("sc_sparse")


#' @export
sc_sparse.default <- function(x, y, out.rm = TRUE, binner = "hex"){
  sc <- scree(x, y, out.rm = out.rm, binner = binner)
  sc_sparse.scree(sc)
}


#' @export
sc_sparse.scree <- function(x, y=NULL, out.rm = FALSE, binner = NULL) {
  #generate vector of MST edges
  mst <- gen_mst(x$del, x$weights)
  sc_sparse.igraph(mst,x)
}



#' @export
sc_sparse.igraph <- function(x, y, out.rm = TRUE, binner = "hex"){
  #input: x is MST, y is scree
  mstmat <- twomstmat(x,y)$lowertri
  edges <- mstmat[which(mstmat>0)]
  #calculate sample size weight
  n = length(x)/500
  w = 0.7 + 0.3/(1+n^2)
  #calculate sparse value
  w*(sort(edges)[floor(0.9*length( edges))])
}


#' Compute skewed scagnostic measure using MST
#'
#' A measure of skewness in the edge lengths of the MST (not in the
#' distribution of the data). It was first defined in Graph Theoretic
#' Scagnostics, Wilkinson, et al. (2005). It is the ratio between
#' the 90th to 50th percentile range and the 10th to 90th percentile range.
#'
#' @inheritParams scree
#' @return A numeric object that gives the plot's skewed score.
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(datasaurus_dozen, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~dataset, ncol=3, scales = "free")
#'   sc_skewed(datasaurus_dozen_wide$dots_x, datasaurus_dozen_wide$dots_y)
#'   sc_skewed(datasaurus_dozen_wide$h_lines_x, datasaurus_dozen_wide$h_lines_y)
#'   sc_skewed(datasaurus_dozen_wide$x_shape_x, datasaurus_dozen_wide$x_shape_y)
#'
#' @export
sc_skewed <- function(x, y, out.rm = TRUE, binner = "hex") UseMethod("sc_skewed")


#' @export
sc_skewed.default <- function(x, y, out.rm = TRUE, binner = "hex"){
  sc <- scree(x, y, out.rm = out.rm, binner = binner)
  sc_skewed.scree(sc)
}


#' @export
sc_skewed.scree <- function(x, y=NULL, out.rm = FALSE, binner = NULL) {
  #generate vector of MST edges
  mst <- gen_mst(x$del, x$weights)
  sc_skewed.igraph(mst, x)
}


#' @export
sc_skewed.igraph <- function(x, y, out.rm = TRUE, binner = "hex"){
  mstmat <- twomstmat(x,y)$lowertri
  edges <- mstmat[which(mstmat>0)]

  # find quantiles
  q10 <- sort(edges)[floor(0.1*length( edges))]
  q50 <- sort(edges)[floor(0.5*length( edges))]
  q90 <- sort(edges)[floor(0.9*length( edges))]

  #calculate sample size weight
  n = length(x)/500
  w = 0.7 + 0.3/(1+n^2)

  # calculate skewed value
  w*((q90-q50)/(q90-q10))

}


#' Compute outlying scagnostic measure using MST
#'
#' A measure of proportion and severity of outliers in the dataset. It was
#' first defined in Graph Theoretic Scagnostics, Wilkinson, et al. (2005).
#' It is calculated by comparing the edge lengths of the outlying points in
#' the MST with the total length of all the edges in the MST.
#'
#' @inheritParams scree
#' @return A numeric object that gives the plot's outlying score.
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   ggplot(datasaurus_dozen, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~dataset, ncol=3, scales = "free")
#'   sc_outlying(datasaurus_dozen_wide$dino_x, datasaurus_dozen_wide$dino_y)
#'   sc_outlying(datasaurus_dozen_wide$dots_x, datasaurus_dozen_wide$dots_y)
#'   sc_outlying(datasaurus_dozen_wide$h_lines_x, datasaurus_dozen_wide$h_lines_y)
#'
#' @export
sc_outlying <- function(x, y, out.rm = TRUE, binner = "hex") UseMethod("sc_outlying")


#' @export
sc_outlying.default <- function(x, y, out.rm = TRUE, binner = "hex"){
  sc <- scree(x, y, out.rm = out.rm, binner = binner)
  sc_outlying.scree(sc)
}


#' @export
sc_outlying.scree <- function(x, y=NULL, out.rm = FALSE, binner = NULL) {
  #generate vector of MST edges
  mst <- gen_mst(x$del, x$weights)
  sc_outlying.igraph(mst, x)
}


#' @export
sc_outlying.igraph <- function(x, y, out.rm = TRUE, binner = "hex"){
  #input: x orig mst (mymst) and y scree object
  #output: outlying mst value

  #make into matrix
  twomst <- twomstmat(x,y)
  mstmat <- twomst$mat
  mstlowertri <- twomst$lowertri

  #identify outliers
  outliers <- outlying_identify(x, y)

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


