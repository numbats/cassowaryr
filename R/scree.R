#' Pre-processing to generate scagnostic measures
#'
#' @param x,y numeric vectors
#' @param binner an optional function that bins the x and y vectors prior
#' to triangulation
#' @param ...  other args
#'
#'
#' @return An object of class "scree" that consists of four elements:
#'  - `del`: the Delauney-Voronoi tesselation from [alphahull::delvor()]
#'  - `weights`: the lengths of each edge in the Delauney triangulation
#'  - `outlier_rm`: cutoff threshold for MST-based outlier detection
#'  - `alpha_param`: the radius or `alpha_param` value calculated from
#'  90th percentile of MST edge lengths
#'
#' @examples
#'
#' x <- runif(100)
#' y <- runif(100)
#' scree(x,y)
#'
#' @importFrom stats median quantile
#' @export
scree <- function(x, y, binner = NULL, alpha_param = NULL, ...) {
  # checks on x,y
  stopifnot(
    is.numeric(x), is.numeric(y), length(x) == length(y)
  )
  # Check if data is a straight line
  if (any(abs(stats::cor(x,y))>1-1*10^-15, !stats::sd(x)>0, !stats::sd(y)>0))
    stop("Data is a perfectly straight line and cannot be analysed")

  # Binner must be a function
  if (!(is.null(binner) | is.function(binner)))
    stop("binner must be a function")

  # cast to a matrix
  xy <- cbind(unitize(x), unitize(y))

  # Remove duplicates
  xy <- xy[!duplicated(round(xy, 10)), , drop = FALSE]

  # Binner function
  if (is.function(binner)) {
    xy <- binner(xy)
  }

  # compute delauney triangulation
  del <- alphahull::delvor(xy)

  # edge weights from the triangulation
  weights <- gen_edge_lengths(del)

  # Compute outlier removal threshold (Wilkinson’s boxplot rule)
  outlier_rm <- psi(weights)

  # Compute alpha radius if not provided (Wilkinson’s suggestion: 90th percentile of MST edge lengths)

  # Generate MST and compute its edge weights
  mst <- gen_mst(del, weights)
  mst_weights <- igraph::E(mst)$weight

  if (is.null(alpha_param)) {
    alpha_param <- stats::quantile(mst_weights, 0.9)
  } else {
    alpha_param <- alpha_param
  }

  structure(
    list(
      del = del,
      weights = weights,
      outlier_rm = outlier_rm,
      alpha_param = alpha_param
    ),
    class = "scree"
  )
}

gen_edge_lengths <- function(del) {
  from_cols <- c("x1", "y1")
  to_cols <- c("x2", "y2")
  sqrt(rowSums((del$mesh[, from_cols] - del$mesh[, to_cols])^2))
}

# rescale input to lie in unit interval
unitize <- function(x, na.rm = TRUE) {
  rng <- range(x, na.rm = na.rm)
  (x - rng[1]) / diff(rng)
}

# This is the edge filter from Wilkinson 05
psi <- function(w, q = c(0.25, 0.75)) {
  q <- stats::quantile(w, probs = q)
  unname(q[2] + 1.5 * diff(q))
}

# rescale input to lie in unit interval
unitize <- function(x, na.rm = TRUE) {
  rng <- range(x, na.rm = na.rm)
  (x - rng[1]) / diff(rng)
}

# This is the edge filter from Wilkinson 05
psi <- function(w, q = c(0.25, 0.75)) {
  q <- stats::quantile(w, probs = q)
  unname(q[2] + 1.5 * diff(q))
}
