####################################################################################
# CLUMPY2: function that gets medians from clusters after
get_median_k <- function(mst, j){
  # get edge to delete
  edges <- igraph::E(mst)[[j]]

  # delete edges
  disjoint_mst <- decompose(igraph::delete_edges(mst, edges))
  print("check")
  # extract stats from disjoint graphs
  gs <- graph_stats(disjoint_mst) |>
    data.frame() |>
    stats::na.omit() # remove single nodes (rows with NA median_edge)

  # K = max_edge from MST with min n
  return(gs[ , "med_edge"])
}


# find max gap
arrange_edges <- sort(mst_weights, decreasing = TRUE)
ind <- which.min(diff(arrange_edges))
# get inter-cluster edges
inter_edges <- utils::head(arrange_edges, ind)

# WRITE FUNCTIO THAT DELETES
# ALL OTHER INTER-CLUSTER EDGES THEN
# FINDS THE GRAPH WITH THE EDGE WE ARE CHECKING
# THEN ANOTHER FUNCTION THAT COMPUTES THE RATIO THING ON THAT GRAPH
inter_edge_ind <- which(mst_weights %in% inter_edges)

# get median value for inter cluster edges
within_edges <- get_median_k(mst, inter_edge_ind)




1 - (mean(within_edges)/mean(inter_edges))


###################################################################################


# PLACEHOLDER FOR CHECKS
# Remove missing values
d <- tibble::tibble(x=x, y=y)
d <- d[stats::complete.cases(d),]
if (length(x) > nrow(d))
  message("WARNING: ", length(x)-nrow(d), " observations in have been removed. \n")
x <- d$x
y <- d$y

if(is.null(sm_list)){
  # this is null when one of the variables is constant after outlier removal
  return(
    dplyr::tibble("outlying"=NA,
                  "stringy"=NA,
                  "striated"=NA,
                  "grid" = NA,
                  "clumpy"=NA,
                  "clumpy2"=NA,
                  "sparse"=NA,
                  "skewed"=NA,
                  "convex"=NA,
                  "skinny"=NA,
                  "monotonic"=NA,
                  "splines"=NA,
                  "dcor"=NA,
                  "stripes"=NA)
  )
}

# require(ggplot2)
# require(tidyr)
# require(dplyr)
#
# data(anscombe_tidy)
# ggplot(anscombe_tidy, aes(x=x, y=y)) +
#   geom_point() +
#   facet_wrap(~set, ncol=2, scales = "free")
#
# draw_mst(anscombe$x1, anscombe$y1)
# draw_mst(anscombe$x2, anscombe$y2, out.rm=TRUE)
# draw_mst(anscombe$x3, anscombe$y3, out.rm=TRUE)
# draw_mst(anscombe$x4, anscombe$y4, out.rm=TRUE)

# sc_stringy(anscombe$x1, anscombe$y1)
# sc_stringy(anscombe$x2, anscombe$y2)
# sc_stringy(anscombe$x3, anscombe$y3)
# sc_stringy(anscombe$x4, anscombe$y4)
#
# sc_striated(anscombe$x1, anscombe$y1)
# sc_striated(anscombe$x2, anscombe$y2)
# sc_striated(anscombe$x3, anscombe$y3)
# sc_striated(anscombe$x4, anscombe$y4)

# test_that("multiplication works", {
#   expect_equal(2 * 2, 4)
# })
