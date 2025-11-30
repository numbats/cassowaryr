# Testing iterative outlier removal
test_that("outlier removal behaves differently for clean vs outlier-augmented data", {

  set.seed(1053)
  # Clean pattern
  x1 <- seq(0, 10, by = 0.01)
  y1 <- sin(2 * x1)

  # Same pattern with two extra outliers
  x2 <- c(x1, -5, 15)
  y2 <- c(y1, -1, 0)

  sc1 <- scree(x1, y1, outlier_rm = FALSE)
  sc2 <- scree(x1, y1, outlier_rm = TRUE)

  sc3 <- scree(x2, y2, outlier_rm = FALSE)
  sc4 <- scree(x2, y2, outlier_rm = TRUE)

  # For clean data, outlier removal should not change the result
  expect_true(isTRUE(all.equal(sc1$alpha,   sc2$alpha)))
  expect_true(isTRUE(all.equal(sc1$weights, sc2$weights)))

  # For data with outliers added, outlier removal should change the result
  expect_false(isTRUE(all.equal(sc3$alpha,   sc4$alpha)))
  expect_false(isTRUE(all.equal(sc3$weights, sc4$weights)))
})



test_that("hexagonal binning changes the scree structure compared to raw points", {

  set.seed(1103)
  x <- seq(0, 10, by = 0.01)
  y <- sin(2 * x)

  # Scree without binning (raw points)
  sc_raw <- scree(x, y, binner = NULL)

  # Scree with hexagonal binning
  sc_hex <- scree(x, y, binner = "hex")

  expect_false(isTRUE(all.equal(sc_raw$weights, sc_hex$weights)))
  expect_false(isTRUE(all.equal(sc_raw$alpha,   sc_hex$alpha)))

  #  we expect fewer edges than with raw points
  expect_gt(length(sc_raw$weights), length(sc_hex$weights))
})
