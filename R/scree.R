#' Pre-processing to generate scagnostic measures
#'
#' @param x,y numeric vectors
#' @param binner an optional function that bins the x and y vectors prior
#' to triangulation
#' @param global_rng an optional vector c(min, max) that specifies the global
#' range of the data set if you want all variables one the same scale
#' @param ...  other args
#'
#' @return An object of class "scree" that consists of three elements:
#'  - `del`: the Delauney-Vornoi tesselation from [alphahull::delvor()]
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
scree <- function(x, y, binner = NULL, global_rng = NULL, ...) {
  # checks on x,y
  stopifnot(
    is.numeric(x), is.numeric(y), length(x) == length(y)
  )

  if (!(is.null(binner) | is.function(binner)))
    stop("binner must be a function")

  # cast to a matrix
  if (is.numeric(global_rng)) xy <- cbind(unitize(x, global_rng), unitize(y, global_rng))
  else xy <- cbind(unitize(x), unitize(y))

  # Check for duplicates and remove
  xy <- xy[!duplicated(xy),]

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
unitize <- function(x, global_rng = NULL, na.rm = TRUE) {
  # set range of x
  rng <- range(x, na.rm = na.rm)

  # if a global range has been specified
  if (is.numeric(global_rng)){

    # stopping criteria
    if (!length(global_rng)==2 | global_rng[1] > global_rng[2]){
      stop("The global range must be a vector with two values, c(min,max).")
    }
    if (rng[1]<global_rng[1] | rng[2]>global_rng[2]){
      stop("The range of x is not a subset of the global range")
    }

    # reset x range
    rng <- global_rng
  }

  # scale 0 to 1
  (x - rng[1]) / diff(rng)
}

# This is the edge filter from Wilkinson 05
psi <- function(w, q = c(0.25, 0.75)) {
  q <- stats::quantile(w, probs = q)
  unname(q[2] + 1.5 * diff(q))
}
