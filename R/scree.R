#' Pre-processing to generate scagnostic measures
#'
#' This function performs the pre-processing requires to calculate the
#' scagnostic measures. This includes the binning, outlier removal, and
#' calculation of the alpha value.
#'
#' @param x,y numeric vectors
#' @param binner an optional function that bins the x and y vectors prior
#' to triangulation
#'  Can be:
#'   \itemize{
#'     \item `"hex"` (default): hexagonal binning following the procedure in the
#'     graph-theoretic scagnostics paper (start 40x40, halve
#'              until <= 250 nonempty cells)
#'     \item `NULL`: no binning (use raw points)
#'     \item a function: user-defined binner
#'     }
#' @param alpha character, numeric, or function. Controls the alpha radius.
#'   Valid character values are:
#'   \itemize{
#'     \item "rahman" (default): Rahman's MST-based middle-50% alpha
#'     \item "q90": 90th percentile of MST edge lengths
#'     \item "omega": graph-theoretic scagnostics alpha
#'   Alternatively:
#'     \item a numeric value giving a fixed alpha
#'     \item a function with no arguments that returns a single numeric alpha
#'     }
#' @param out.rm logical; if TRUE, iteratively trim large MST edges, If FALSE
#' the scagnostics will be computed on the entire data set with no outlier
#' removal.
#' @param ...  other args
#'
#' @return An object of class "scree" that consists of three elements:
#'  - `del`: the Delauney-Voronoi tesselation from [alphahull::delvor()]
#'  - `weights`: the lengths of each edge in the Delauney triangulation
#'  - `alpha`: the radius or `alpha` value that will be used to generate the
#'  alphahull
#'
#' @examples
#' set.seed(232)
#'
#' x <- runif(1000)
#' y <- runif(1000)
#'
#' # make scree
#' sc0 <- scree(x,y)
#' sc1 <- scree(x,y, out.rm = FALSE)  # no outlier removal
#' sc2 <- scree(x, y, binner = NULL) #  no hexagonal binning
#'
#' # see the difference made by binning out out.rm
#' draw_mst(sc0)
#' draw_mst(sc1)
#' draw_mst(sc2)
#'
#' @export
scree <- function(x, y, out.rm = TRUE, binner = "hex",
                  alpha = "rahman", ...) {
  # Check if data is numeric and vectors are equal length
  stopifnot(
    is.numeric(x), is.numeric(y), length(x) == length(y)
  )
  # cast to a matrix
  xy <- cbind(unitize(x), unitize(y))

  # Check for duplicates using same method as interp (otherwise get error)
  xy[!is_dupe(x,y),]

  # Bin data
  xy <- get_binned_matrix(xy, binner)

  # Check for perfect colinearity
  if (any(abs(stats::cor(x, y)) > 1-1*10^-15, # straight line
          !stats::sd(x)>0, # x is constant
          !stats::sd(y)>0) # y is constant
  ){
    message("Data is a perfectly straight line and cannot be analysed. ",
            "The Delaunay triangulation will be computed with a small amount of jitter.")
    # make matrix of very small noise
    noise_matrix <- matrix(rnorm(prod(dim(xy)), 0, 1*10^-5),
                           dim(xy)[1], dim(xy)[2])
    # add noise to original matrix
    xy <- xy + noise_matrix
  }

  # Outlier removal & delauney triangulation
  if (out.rm) {
    del <- outlier_removal(xy)
  } else{
    del <- alphahull::delvor(xy)
  }
  # Full graph edge weights
  weights <- gen_edge_lengths(del)
  # Alpha value
  alpha_value <- get_numeric_alpha(alpha, del, weights)


  structure(
    list(
      del = del,
      weights  = weights,
      alpha = alpha_value
    ),
    class = "scree"
  )
}

