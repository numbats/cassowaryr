######################### OUTLIER REMOVAL FUNCTIONS #########################
outlier_removal <- function(){
  # Generate MST and compute its edge weights
  mst <- gen_mst(del, weights)
  mst_weights <- igraph::E(mst)$weight

  # remove outliers iteratively (based on "Scagnostics Distributions" paper)
  if (outlier_rm) {
    repeat {
      w <- psi(mst_weights)
      long_edges <- which(mst_weights > w)

      if (length(long_edges) == 0L) break

      long_vertices <- igraph::ends(mst, long_edges)
      outlier_vertices <- unique(as.integer(long_edges))

      if (length(outlier_vertices) == 0L ||
          nrow(xy) - length(outlier_vertices) < 3L) {
        message("Outlier removal stopped: removing these outliers would leave fewer than 3 points (too few vertices to continue).")
        break
      }

      xy <- xy[-outlier_vertices, , drop = FALSE]

      # recompute del, weights, mst, mst_weights on updated xy
      del     <- alphahull::delvor(xy)
      weights <- gen_edge_lengths(del)
      mst     <- gen_mst(del, weights)
      mst_weights <- igraph::E(mst)$weight
    }
  }
}


outlying_identify <- function(mst, sc){
  #input: takes a mst and scree
  #output: rown number of

  #get matrix and upper triangular matrix
  matlist <- twomstmat(mst,sc)
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

original_and_robust <- function(x, y){
  #input: data for 2 variables x and y
  #output: list of scree and MST objects

  #construct original scree and MST
  sc_original <- scree(x, y)
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
    if(stats::sd(new_x) == 0 | stats::sd(new_y) == 0) return(NULL)
    #recalculate scree and MST
    sc_robust <- scree(new_x, new_y)
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
