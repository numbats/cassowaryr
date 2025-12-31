

get_binned_matrix <- function(xy, binner){
  # Check binner
  if (!is.null(binner) && !is.character(binner) && !is.function(binner)) {
    stop("binner must be NULL, 'hex', or a function")
  }

  # Binning
  if (is.null(binner)) {

    # default: no binning, do nothing

  } else if (is.character(binner)) {

    binner_choice <- match.arg(binner, "hex")

    if (binner_choice == "hex") {
      xy <- hex_binner(xy, ...)
    }

  } else if (is.function(binner)) {

    xy <- binner(xy, ...)

    if (!is.matrix(xy) || ncol(xy) != 2) {
      stop("User-supplied binner must return a two-column matrix of (x, y) points.")
    }
  }
  xy
}

# Hexagonal binning as in the graph-theoretic scagnostics paper
hex_binner <- function(xy, xbins = 40, max_cells = 250) {
  if (!requireNamespace("hexbin", quietly = TRUE)) {
    stop("Package 'hexbin' must be installed to use binner = 'hex'.")
  }

  current_xbins <- xbins
  repeat {
    hb <- hexbin::hexbin(xy[, 1], xy[, 2], xbins = current_xbins)
    n_cells <- length(hb@count)  # number of non-empty hex cells
    print(n_cells)
    if (n_cells <= max_cells || current_xbins <= 1) {
      break
    }

    # reduce bin size by half and rebin
    current_xbins <- max(1L, floor(current_xbins / 2))
  }

  centers <- hexbin::hcell2xy(hb)
  cbind(centers$x, centers$y)
}



