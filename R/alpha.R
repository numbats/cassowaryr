# internal function for scree that returns numeric alpha value
get_numeric_alpha <- function(alpha, del, weights){
  if (is.character(alpha)) {
    # Most definitions depend on the MST edge lengths
    mst <- gen_mst(del, weights)
    mst_weights <- igraph::E(mst)$weight

    # get required params
    alpha_choice <- match.arg(alpha, c( "rahman", "q90", "omega"))
    n <- nrow(xy)

    # calculate set alpha value
    alpha_value <- switch(
      alpha_choice,
      rahman = alpha_rahman(mst_weights),
      q90    = alpha_q90(mst_weights),
      omega  = alpha_omega(mst_weights)
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
  return(alpha_value)
}

# Alpha value using 90th percentile of MST edge length
alpha_q90 <- function(mst_weights) {
  unname(stats::quantile(mst_weights, 0.9))
}

# Rahman's suggested MST-based alpha
alpha_rahman <- function(mst_weights) {
  n <- length(mst_weights) + 1
  q <- stats::quantile(mst_weights, probs = c(0.25, 0.75))
  middle_edges <- mst_weights[mst_weights >= q[1] & mst_weights <= q[2]]
  sqrt(sum(middle_edges) / n)
}

# Alpha value suggested in "Graph Theoretic Scagnostics" paper
alpha_omega <- function(mst_weights) {
  psi(mst_weights)
}
