# make data function
set.seed(232)

x <- runif(1000)
y <- runif(1000)

# save correct outputs
# sc0 <- scree(x,y, binner = NULL, out.rm = FALSE)
# sc1 <- scree(x,y, binner = NULL)  # no binning, remove outliers
# sc2 <- scree(x, y, out.rm = FALSE) #  hexagonal binning
# sc3 <- scree(x, y) # both
# saveRDS(sc0, file = "tests/testthat/objects/sc0.rds")
# saveRDS(sc1, file = "tests/testthat/objects/sc1.rds")
# saveRDS(sc2, file = "tests/testthat/objects/sc2.rds")
# saveRDS(sc3, file = "tests/testthat/objects/sc3.rds")

# read in correct outputs
sc0 <- readRDS("objects/sc0.rds")
sc1 <- readRDS("objects/sc1.rds")
sc2 <- readRDS("objects/sc2.rds")
sc3 <- readRDS("objects/sc3.rds")

# Testing iterative outlier removal
test_that("test scree function works", {
  # just check functions have expected output
  expect_equal(scree(x,y, binner = NULL, out.rm = FALSE), sc0)
  expect_equal(scree(x,y, binner = NULL), sc1)
  expect_equal(scree(x, y, out.rm = FALSE), sc2)
  expect_equal(scree(x, y), sc3)

  #  we expect fewer edges than with raw points
  expect_gt(length(sc0$weights), length(sc1$weights))
  expect_gt(length(sc0$weights), length(sc2$weights))
})
