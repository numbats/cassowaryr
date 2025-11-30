# Load Libraries
library(dplyr)

# Set Seed
set.seed(4213)

# Generate circle of data for test
data <- tibble(r = runif(20000, 0,1),
               theta = runif(20000, 0, 2*pi)) %>%
  mutate(x = r*cos(theta),
         y = r*sin(theta))

#plot(data$x,data$y)

# Test Skinny Scagnostic
test_that("Skinny Scagnostic", {
  expect_lt(sc_skinny(data$x,data$y), 0.1)
})

# Test Convex Scagnostic
test_that("Convex Scagnostic", {
  expect_gt(sc_convex(data$x,data$y), 0.95)
})


# additional tests for skinny index
test_that("Skinny scagnostic", {
  set.seed(1200)

  x  <- seq(0, 1, length.out = 50)

  y1 <- x
  y2 <- 2 * x

  s1 <- suppressWarnings(sc_skinny(x, y1))
  s2 <- suppressWarnings(sc_skinny(x, y2))

  # we expect exactly 1
  expect_equal(s1, 1, tolerance = 1e-8)
  expect_equal(s2, 1, tolerance = 1e-8)
})


