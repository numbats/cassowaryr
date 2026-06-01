test_that("Test diagnose_outliers", {

  set.seed(212)

  data <- tibble::tibble(
    id = seq_len(35),
    group = rep(c("A", "B"), length.out = 35),
    var1 = c(rnorm(30), 8, 9, -8, -9, 10),
    var2 = c(rnorm(30), 8, -9, 8, -9, 10)
  )

  res <- diagnose_outliers(data, var1, var2)

  expect_equal(nrow(res), nrow(data))
  expect_true("outlier_status" %in% names(res))
  expect_true(all(res$outlier_status %in% c("Kept", "Removed")))
})
