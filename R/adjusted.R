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
#' require(ggplot2)
#' require(dplyr)
#'
#' # plot features
#' ggplot(features, aes(x=x, y=y)) +
#'    geom_point() +
#'    facet_wrap(~feature, ncol = 5, scales = "free")
#'
#' # calculate clumpy2 on all features
#' features |>
#'   group_by(feature) |>
#'   summarise(clumpy2 = sc_clumpy2(x,y))
#'
#' sc_clumpy2(datasaurus_dozen_wide$dots_x, datasaurus_dozen_wide$dots_y)
#'
#' data <- features |> filter(feature == "clusters")
#' x <- data$x
#' y <- data$y
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
sc_clumpy2.scree <- function(x, y=NULL, out.rm = FALSE, binner =  "hex") {
  mst <- gen_mst(x$del, x$weights)
  sc_clumpy2.igraph(mst,x)
}


#' @export
sc_clumpy2.igraph <- function(x, y=NULL, out.rm = TRUE, binner =  "hex"){
  # rename x to mst
  mst <- x
  #set stringy penalty
  vertex_counts <- igraph::degree(x)
  stringy <- sum(vertex_counts == 2) / (length(vertex_counts) - sum(vertex_counts == 1))
  stringy_pen <- ifelse(stringy>0.9, (1-stringy), 1)

  # find max gap
  mst_weights <- igraph::E(mst)$weight
  arrange_edges <- sort(mst_weights, decreasing = TRUE)
  ind <- which.min(diff(arrange_edges))

  # get inter-cluster edges
  inter_edges <- utils::head(arrange_edges, ind)
  inter_edge_ind <- which(mst_weights %in% inter_edges)
  # add edge attribute
  mst <- igraph::set_edge_attr(mst, "inter", index = inter_edge_ind, TRUE)

  # calculate clumpy
  clumpy_vector <- t(sapply(seq(length(inter_edge_ind)),
                          # Function: given edge index, find clumpy value
                          function(x) {
                            # get stats from smaller of two clusters
                            cluster_graph <- connected_subgraph(mst, inter_edge_ind, x)
                            ntotal <- igraph::vcount(cluster_graph)

                            # two clusters cannot have 3 points or fewer
                            if(ntotal<=3){
                              final_clumpy = 0
                            } else{
                              # get index of edge of interest (changed with subgraph)
                              e_ind <- which(igraph::E(cluster_graph)$inter)
                              graph_stats <- get_graph_feature(cluster_graph, e_ind, c("max_edge", "n"))
                              # Calculate clumpy metric for this value of j
                              clumpy <- 1 - (graph_stats$max_edge/mst_weights[inter_edge_ind[x]])

                              # calculate uneven cluster pentalty
                              nsmall <- graph_stats$n
                              size_penalty <- sqrt(2*nsmall/ntotal)

                              # Final clumpy measure for j (including weighting)
                              final_clumpy <- ntotal*clumpy*size_penalty
                            }

                            # Take number of points so that we can weight final value
                            c(clumpy = final_clumpy, n = ntotal)

                          }
  ))
  # return clump value and finish weighting
  clumpy_vector <- clumpy_vector |> data.frame()
  stringy_pen*(sum(clumpy_vector$clumpy)/sum(clumpy_vector$n))
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

# get mst with all edges deleted
# input: mst, full edge index, index of edge j
# output: subgraph of two clusters that are connected by edge j
connected_subgraph <- function(mst, inter_edge_ind, j){
  # get fully disconnected mst with all edges deleted
  edges_all <- igraph::E(mst)[[inter_edge_ind]]
  disjoint_all <- igraph::delete_edges(mst, edges_all)
  group_all <- igraph::groups(igraph::components(disjoint_all))

  # get mst with all but index j edge deleted
  e_ind <- inter_edge_ind[-j] # remove other edges
  edges_but_one <- igraph::E(mst)[[e_ind]]
  disjoint_but_one <-igraph::delete_edges(mst, edges_but_one)
  group_but_one <- igraph::groups(igraph::components(disjoint_but_one))

  # get vertex indicies for mismatched vertice list
  verts <- group_but_one[which(!(group_but_one %in%  group_all))]
  verts <- unname(unlist(verts))
  # make a subgraph with those verticies
  igraph::subgraph(disjoint_but_one, verts)
}


