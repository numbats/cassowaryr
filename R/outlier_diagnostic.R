#' Diagnose outlier removal for one variable pair
#'
#' Identifies which observations are kept or removed by the outlier removal
#' step for a chosen pair of numeric variables.
#'
#' @param data A data frame or tibble.
#' @param x,y Two numeric columns to use for the outlier removal diagnostic.
#'
#' @return A data frame with the same rows and columns as `data`, plus an
#'   `outlier_status` column. Values are `"Kept"` or `"Removed"`.
#'
#' @examples
#' data <- data.frame(
#'   id = 1:10,
#'   var1 = c(1, 2, 3, 4, 5, 6, 7, 20, 21, 22),
#'   var2 = c(1, 2, 3, 4, 5, 6, 7, 20, -20, 22)
#' )
#'
#' diagnose_outliers(data, var1, var2)
#'
#' @export
diagnose_outliers <- function(data, x, y) {

  x_vals <- dplyr::pull(data, {{ x }})
  y_vals <- dplyr::pull(data, {{ y }})


  if (!is.numeric(x_vals)) {
    stop("`x` must be numeric.")
  }

  if (!is.numeric(y_vals)) {
    stop("`y` must be numeric.")
  }

  x_unit <- unitize(x_vals)
  y_unit <- unitize(y_vals)

  keep_idx <- is_dupe(x_vals, y_vals)

  xy_mat <- cbind(x_unit, y_unit)
  xy_mat <- xy_mat[keep_idx, , drop = FALSE]

  final_del <- outlier_removal(xy_mat)


  kept_tbl <- tibble::tibble(
    x_unit = final_del$x[, 1],
    y_unit = final_del$x[, 2],
    outlier_status = "Kept"
  )


  data |>
    dplyr::mutate(
      x_unit = x_unit,
      y_unit = y_unit
    ) |>
    dplyr::left_join(kept_tbl, by = c("x_unit", "y_unit")) |>
    dplyr::mutate(
      outlier_status = dplyr::if_else(
        keep_idx & !is.na(.data$outlier_status), "Kept", "Removed")
    ) |>
    dplyr::select(-x_unit, -y_unit)
}
