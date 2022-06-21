#' Pre-processing to generate scagnostic measures
#'
#' @param x,y numeric vectors
#' @param binner an optional function that bins the x and y vectors prior
#' to triangulation
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
scree <- function(x, y, binner = NULL, ...) {
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
  dupes <- paste(x, y, sep =",")
  xy <- xy[!duplicated(dupes),]

  # Binner function
  if (is.function(binner)) {
    xy <- binner(xy)
  }

  # compute delauney triangulation
  del <- alphahull::delvor(xy)

  # edge weights from the triangulation
  weights <- gen_edge_lengths(del)

  # alpha estimator
  alpha <- psi(weights)

  structure(
    list(
      del = del,
      weights  = weights,
      alpha = alpha
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
