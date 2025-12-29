
######################### UTILITY FUNCTIONS ##############################

# remove duplicates that would trigger an error from the interp package
is_dupe <- function(x, y){
  xrnd <- round(unitize(x), digits = 10)
  yrnd <- round(unitize(y), digits = 10)
  xy_chars <- paste(xrnd, yrnd, sep =",")
  !duplicated(xy_chars)
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

