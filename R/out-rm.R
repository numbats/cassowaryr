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

    # Finally, we need to check if outlier removal created straight lines
    corr <- suppressWarnings(abs(stats::cor(xy[,1], xy[,2])) > 1-1*10^-15)
    # In which case we need to add noise to the xy matrix
    if (any(corr, # straight line
            !stats::sd(xy[,1])>0, # x is constant
            !stats::sd(xy[,2])>0) # y is constant
    ){
      message("Outlier removal produced data that is in a perfectly straight line. ",
              "The Delaunay triangulation will be computed with a small amount of jitter.")
      # make matrix of very small noise
      noise_matrix <- matrix(rnorm(prod(dim(xy)), 0, 1*10^-5),
                             dim(xy)[1], dim(xy)[2])
      # add noise to original matrix
      xy <- xy + noise_matrix
    }
  }
}

