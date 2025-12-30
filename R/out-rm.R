outlying_identify <- function(mst, mst_weights){
  #input: takes a mst and scree
  #output: row index of outliers

  # get number of observations
  n <- length(mst_weights) + 1

  # get edge matrix of MST
  mstmat <- matrix(mst[], nrow = n)

  #calculate omega value
  w <- psi(mst_weights)

  # set values above omega to 0 in matrix to find outlying
  mstmat_check <- mstmat
  mstmat_check[mstmat>w]=0

  #row sum of matrix, if 0 all edges are 0 then it is an outlier
  rowsum <- mstmat_check%*%rep(1, length(mstmat_check[1,]))

  if(sum(rowsum==0) == 0){
    # if no outliers return NULL
    return(NULL)
  } else{
    # Otherwise return the index of the outlying observation
    return(which(rowsum==0))
  }
}

outlier_removal <- function(xy){
  # remove outliers iteratively (based on "Scagnostics Distributions" paper)
  repeat {
    # Get del, MST and weights
    del <- alphahull::delvor(xy)
    weights <- gen_edge_lengths(del)
    mst <- gen_mst(del, weights)
    mst_weights <- igraph::E(mst)$weight

    # get outliers
    outliers <- outlying_identify(mst, mst_weights)
    # if no outliers, just return current del
    if(is.null(outliers)){
      return(del)
    } else{
      # otherwise we remove outliers and recompute
      xy <- xy[-outliers,]
    }
  }
}

