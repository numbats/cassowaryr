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
