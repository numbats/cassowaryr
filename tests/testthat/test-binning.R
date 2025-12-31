set.seed(1989) # good album
x <- rnorm(1000)
y <- rnorm(1000)

# get matrix object
xy <- cbind(unitize(x), unitize(y))

# get binned object
# code that made hexbin_xy
# hexbin_xy <- get_binned_matrix(xy, "hex")
# saveRDS(hexbin_xy,
#         file = "tests/testthat/objects/hexbin_xy.rds")
hexbin_xy <- readRDS("tests/testthat/objects/hexbin_xy.rds")


test_that("Test get_binned_matrix", {
  # no binning
  expect_equal(get_binned_matrix(xy, NULL), xy)
  # hex binning from scagnostics package
  expect_equal(get_binned_matrix(xy, "hex"), hexbin_xy)
})

