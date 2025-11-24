#' Pre-processing to generate scagnostic measures
#'
#' @param x,y numeric vectors
#' @param binner an optional function that bins the x and y vectors prior
#' to triangulation
#' @param alpha character, numeric, or function. Controls the alpha radius.
#'   Valid character values are:
#'     - "rahman" (default): Rahman's MST-based middle-50% alpha
#'     - "q90": 90th percentile of MST edge lengths
#'     - "omega": graph-theoretic scagnostics alpha
#'   Alternatively:
#'     - a numeric value giving a fixed alpha
#'     - a function with no arguments that returns a single numeric alpha
#' @param ...  other args
#'
#' @return An object of class "scree" that consists of three elements:
#'  - `del`: the Delauney-Voronoi tesselation from [alphahull::delvor()]
#'  - `weights`: the lengths of each edge in the Delauney triangulation
#'  - `alpha`: the radius or `alpha` value that will be used to generate the
#'  alphahull
#'
#' @examples
#'
#' x <- runif(100)
#' y <- runif(100)
#' scree(x,y)
#'
#' @export
scree <- function(x, y, binner = NULL, alpha = c("rahman", "q90", "omega"), ...) {
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

  # Check for duplicates and remove
  # (had to cut off at 15 digits otherwise shull spits error)
  xrnd <- round(unitize(x), digits = 10)
  yrnd <- round(unitize(y), digits = 10)
  dupes <- paste(xrnd, yrnd, sep =",")
  xy <- xy[!duplicated(dupes),]

  # Binner function
  if (is.function(binner)) {
    xy <- binner(xy)
  }

  # compute delauney triangulation
  del <- alphahull::delvor(xy)

  # edge weights from the triangulation
  weights <- gen_edge_lengths(del)

  # Generate MST and compute its edge weights
  mst <- gen_mst(del, weights)
  mst_weights <- igraph::E(mst)$weight

  n <- nrow(xy)


  if (is.character(alpha)) {

    alpha_choice <- match.arg(alpha, c("rahman", "q90", "omega"))

    alpha_value <- switch(
      alpha_choice,
      rahman = alpha_rahman(mst_weights, n),
      q90    = alpha_q90(mst_weights),
      omega  = alpha_omega(weights)
    )

  } else if (is.numeric(alpha)) {

    # user-specified numeric value
    alpha_value <- alpha

  } else if (is.function(alpha)) {

    # user-supplied alpha estimator
    alpha_value <- alpha()

    if (!is.numeric(alpha_value) || length(alpha_value) != 1) {
      stop("User-supplied alpha function must return a single numeric value.")
    }

  } else {
    stop("alpha must be a character string, numeric value, or function.")
  }


  structure(
    list(
      del = del,
      weights  = weights,
      alpha = alpha_value
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

# Alpha value using 90th percentile of MST edge length
alpha_q90 <- function(mst_weights) {
  stats::quantile(mst_weights, 0.9)
}

# Rahman's suggested MST-based alpha
alpha_rahman <- function(mst_weights, n) {
  q <- stats::quantile(mst_weights, probs = c(0.25, 0.75))
  middle_edges <- mst_weights[mst_weights >= q[1] & mst_weights <= q[2]]
  sqrt(sum(middle_edges) / n)
}

# Alpha value suggested in "Graph Theoretic Scagnostics" paper
alpha_omega <- function(weights) {
  psi(weights)
}

