
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
#' require(ggplot2)
#' require(dplyr)
#'
#' # plot the features data
#' ggplot(features, aes(x=x, y=y)) +
#'   geom_point() +
#'   facet_wrap(~feature, ncol = 5, scales = "free")
#'
#' # calculate using tidy code
#' features |>
#'  group_by(feature) |>
#'  summarise(stringy = sc_stringy(x,y))
#'
#' # using just vectors of points
#' x <- datasaurus_dozen_wide$star_x
#' y <- datasaurus_dozen_wide$star_y
#'
#' # plot it
#' ggplot() +
#'   geom_point(aes(x = x, y = y))
#'
#' # calculate using vectors
#' sc_stringy(x, y)
#' @export
sc_stringy <- function(x, y, out.rm = TRUE, binner = "hex") UseMethod("sc_stringy")


#' @export
sc_stringy.default <- function(x, y, out.rm = TRUE, binner = "hex"){
  #input: x and y are vectors
  sc <- scree(x, y, out.rm = out.rm, binner = binner)

  sc_stringy.scree(sc)
}


#' @export
sc_stringy.scree <- function(x, y=NULL, out.rm = TRUE, binner = "hex") {
  #input: x is a scree, no y
  mst <- gen_mst(x$del, x$weights)
  sc_stringy.igraph(mst)
}


#' @export
sc_stringy.igraph <- function(x, y=NULL, out.rm = TRUE, binner = "hex") {
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
#' require(ggplot2)
#' require(dplyr)
#'
#' # plot the features data
#' ggplot(features, aes(x=x, y=y)) +
#'   geom_point() +
#'   facet_wrap(~feature, ncol = 5, scales = "free")
#'
#' # calculate using tidy code
#' features |>
#'  group_by(feature) |>
#'  summarise(striated = sc_striated(x,y)) |>
#'  arrange(striated)
#'
#' # using just vectors of points
#' x <- datasaurus_dozen_wide$v_lines_x
#' y <- datasaurus_dozen_wide$v_lines_y
#'
#' # plot it
#' ggplot() +
#'   geom_point(aes(x = x, y = y))
#'
#' # calculate scagnostic
#' sc_striated(x, y)
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
  # only MST measure that still takes the scree, will figure out
  # code that doesn't need it later
  sc_striated.igraph(mst, x)
}



#' @export
sc_striated.igraph <- function(x, y, out.rm = TRUE, binner = "hex"){
  # input: x is MST, y is scree
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
#' require(ggplot2)
#' require(dplyr)
#'
#' # plot the feature
#' ggplot(features, aes(x=x, y=y)) +
#'    geom_point() +
#'    facet_wrap(~feature, ncol = 5, scales = "free")
#'
#' # calculate using tidy code
#' features |>
#'   group_by(feature) |>
#'   summarise(clumpy = sc_clumpy(x,y))
#'
#' # using two vectors
#' x <- datasaurus_dozen_wide$slant_up_x
#' y <- datasaurus_dozen_wide$slant_up_y
#'
#' # plot it
#' ggplot() +
#'   geom_point(aes(x = x, y = y))
#'
#' # calculate using vectors
#' sc_clumpy(x, y)
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

  sc_clumpy.igraph(mst)
}



#' @export
sc_clumpy.igraph <- function(x, y=NULL, out.rm = TRUE, binner = "hex"){
  # Rename x to mst for code clarity
  mst <- x

  # get all edge weights
  mst_weights <- igraph::E(mst)$weight

  # set up clumpy vector
  n <- length(mst_weights)

  # for each j, calculate clumpy
  clumpy <- sapply(seq(n),
         # Function: given edge index, find clumpy value
         function(x) {
           length_j <- mst_weights[x]
           max_length_k <- get_graph_feature(mst, x, "max_edge")
           1 - (max_length_k/length_j)
         }
  )
  # return clumpy measure
  max(clumpy)
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
#' require(ggplot2)
#' require(dplyr)
#'
#' # plot the feature
#' ggplot(features, aes(x=x, y=y)) +
#'    geom_point() +
#'    facet_wrap(~feature, ncol = 5, scales = "free")
#'
#' # calculate using tidy code
#' features |>
#'   group_by(feature) |>
#'   summarise(sparse = sc_sparse(x,y))
#'
#' # using two vectors
#' x <- datasaurus_dozen_wide$dots_x
#' y <- datasaurus_dozen_wide$dots_y
#'
#' # plot it
#' ggplot() +
#'   geom_point(aes(x = x, y = y))
#'
#' # calculate using vectors
#' sc_sparse(x, y)
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
  # pass to MST version of function
  sc_sparse.igraph(mst)
}



#' @export
sc_sparse.igraph <- function(x, y=NULL, out.rm = TRUE, binner = "hex"){
  # get mst edge weights
  mst_weights <- igraph::E(x)$weight

  #sparse value
  sort(mst_weights)[floor(0.9*length(mst_weights))]
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
#' require(ggplot2)
#' require(dplyr)
#'
#' # plot the feature
#' ggplot(features, aes(x=x, y=y)) +
#'    geom_point() +
#'    facet_wrap(~feature, ncol = 5, scales = "free")
#'
#' # calculate using tidy code
#' features |>
#'   group_by(feature) |>
#'   summarise(skewed = sc_skewed(x,y))
#'
#' # using two vectors
#' x <- datasaurus_dozen_wide$away_x
#' y <- datasaurus_dozen_wide$away_y
#'
#' # plot it
#' ggplot() +
#'   geom_point(aes(x = x, y = y))
#'
#' # calculate using vectors
#' sc_skewed(x, y)
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

  # pass to MST version of function
  sc_skewed.igraph(mst)
}


#' @export
sc_skewed.igraph <- function(x, y=NULL, out.rm = TRUE, binner = "hex"){
  # rename for clarity
  mst_weights <- igraph::E(x)$weight

  # find quantiles
  q10 <- sort(mst_weights)[floor(0.1*length(mst_weights))]
  q50 <- sort(mst_weights)[floor(0.5*length(mst_weights))]
  q90 <- sort(mst_weights)[floor(0.9*length(mst_weights))]

  # Skewed calc
  (q90-q50)/(q90-q10)

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
#' require(ggplot2)
#' require(dplyr)
#'
#' # plot the feature
#' ggplot(features, aes(x=x, y=y)) +
#'    geom_point() +
#'    facet_wrap(~feature, ncol = 5, scales = "free")
#'
#' # calculate using tidy code
#' features |>
#'   group_by(feature) |>
#'   summarise(outlying = sc_outlying(x,y))
#'
#' # using two vectors
#' x <- datasaurus_dozen_wide$away_x
#' y <- datasaurus_dozen_wide$away_y
#'
#' # plot it
#' ggplot() +
#'   geom_point(aes(x = x, y = y))
#'
#' # calculate scag
#' sc_outlying(x, y)
#' @export
sc_outlying <- function(x, y, binner = "hex") UseMethod("sc_outlying")


#' @export
sc_outlying.default <- function(x, y, binner = "hex"){
  sc <- scree(x, y, out.rm = FALSE, binner = binner)
  sc_outlying.scree(sc)
}


#' @export
sc_outlying.scree <- function(x, y = NULL, binner = "hex") {
  #generate vector of MST edges and weights
  mst <- gen_mst(x$del, x$weights)

  # pass to MST version of function
  sc_outlying.igraph(mst)
}


#' @export
sc_outlying.igraph <- function(x, y = NULL, binner = "hex"){
  # rename input (x = mst) for code clarity
  mst <- x

  # get mst edges
  mst_weights <- igraph::E(mst)$weight

  # identify outliers
  outliers <- outlying_identify(mst, mst_weights)

  # get mst as edge matrix
  n <- length(mst_weights) + 1
  mst_matrix <- matrix(mst[], nrow = n)

  # calculate outlying value
  # get sum of edges of outlying points
  outlying_edges <- sum(mst_matrix[outliers,])

  # Get double count from outliers connected to outliers
  double_count <- sum(mst_matrix[as.matrix(expand.grid(outliers,outliers))])

  # subtract edges that were counted twice
  sum_outlying <- outlying_edges - 0.5*double_count

  # get the sum total of all edge lengths in MST (*0.5 for matrix symmtetry)
  sum_edges <-  0.5*sum(mst_matrix)

  # final outlying value
  sum_outlying/sum_edges
}

# function that generates MST from delauney triangulation
gen_mst <- function(del, weights) {
  edges <- del$mesh[, c("ind1", "ind2")]
  graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
  graph <- igraph::set_edge_attr(graph, "weight", value = weights)
  igraph::mst(graph, weights =  igraph::E(graph)$weight)
}


# CLUMPY:
# function that get the feature from the smaller of two clusters
# after an edge has been deleted
get_graph_feature <- function(mst, j, feature){
  # get edge to delete
  edge <- igraph::E(mst)[[j]]

  # delete edge
  disjoint_mst <- decompose(igraph::delete_edges(mst, edge))

  # extract stats from disjoint graphs
  gs <- graph_stats(disjoint_mst) |>
    data.frame() |>
    stats::na.omit() # remove single nodes (rows with NA median_edge)

  # K = max_edge from MST with min n
  return(gs[which.min(gs$n), feature])
}


# CLUMPY: function that extracts statistics from disjoint graphs
graph_stats <- function(dg){
  suppressWarnings( #suppress warning for edges by themselves
    t(sapply(seq_along(dg),
             function(x) c(index = x,
                           n = vcount(dg[[x]]),
                           max_edge = max(igraph::E(dg[[x]])$weight),
                           med_edge = unname(quantile(igraph::E(dg[[x]])$weight,0.5))
             )
    )
    )
  )
}


