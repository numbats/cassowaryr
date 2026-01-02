# make data function
set.seed(232)

x <- runif(1000)
y <- runif(1000)

# save correct outputs
# sc0 <- scree(x,y)
# sc1 <- scree(x,y, out.rm = TRUE)  # no binning, remove outliers
# sc2 <- scree(x, y, binner = "hex") #  hexagonal binning
# sc3 <- scree(x, y, binner = "hex", out.rm = TRUE) # both
# saveRDS(sc0, file = "tests/testthat/objects/sc0.rds")
# saveRDS(sc1, file = "tests/testthat/objects/sc1.rds")
# saveRDS(sc2, file = "tests/testthat/objects/sc2.rds")
# saveRDS(sc3, file = "tests/testthat/objects/sc3.rds")

# read in correct outputs
sc0 <- readRDS("tests/testthat/objects/sc0.rds")
sc1 <- readRDS("tests/testthat/objects/sc1.rds")
sc2 <- readRDS("tests/testthat/objects/sc2.rds")
sc3 <- readRDS("tests/testthat/objects/sc3.rds")

# Testing iterative outlier removal
test_that("test scree function works", {
  # just check functions have expected output
  expect_equal(scree(x,y), sc0)
  expect_equal(scree(x,y, out.rm = TRUE), sc1)
  expect_equal(scree(x, y, binner = "hex"), sc2)
  expect_equal(scree(x, y, binner = "hex", out.rm = TRUE), sc3)

  #  we expect fewer edges than with raw points
  expect_gt(length(sc0$weights), length(sc1$weights))
  expect_gt(length(sc0$weights), length(sc2$weights))
})
