#' Pre-processing to generate scagnostic measures
#'
#' @param x,y numeric vectors
#' @param binner an optional function that bins the x and y vectors prior
#' to triangulation
#'  Can be:
#'   \itemize{
#'     \item `NULL`: no binning (use raw points)
#'     \item `"hex"` : hexagonal binning following the procedure in the
#'     graph-theoretic scagnostics paper (start 40x40, halve
#'              until <= 250 nonempty cells)
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
#' @param outlier_rm logical; if TRUE, iteratively trim large MST edges
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
#' sc1 <- scree(x,y, outlier_rm = TRUE)  # no binning, remove outliers
#' sc2 <- scree(x, y, binner = "hex") #  hexagonal binning
#'
#' # see the difference made by binning out out.rm
#' draw_mst(sc0)
#' draw_mst(sc1)
#' draw_mst(sc2)
#'
#' @export
scree <- function(x, y, outlier_rm = FALSE, binner = NULL,
                  alpha = "rahman", ...) {
  # CHECKS
  # Numeric and equal length
  stopifnot(
    is.numeric(x), is.numeric(y), length(x) == length(y)
  )

  # Data is a straight line
  if (any(abs(stats::cor(x,y))>1-1*10^-15, !stats::sd(x)>0, !stats::sd(y)>0)){
  warning("Data is a perfectly straight line and cannot be analysed; ",
          "Delaunay triangulation is not computed.")
    sc <- list(
      del     = NULL,
      weights = NULL,
      alpha   = NA_real_
    )
    class(sc) <- "scree"
    return(sc)
  }

  # cast to a matrix
  xy <- cbind(unitize(x), unitize(y))

  # Check for duplicates using same method as interp (otherwise get error)
  xy[!is_dupe(x,y),]

  # Bin data
  xy <- get_binned_matrix(xy, binner)

  # Outlier removal & delauney triangulation
  if (outlier_rm) {
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




