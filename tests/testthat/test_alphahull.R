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
test_that("Stringy Scagnotist", {
  expect_lt(sc_skinny(data$x,data$y), 0.1)
})

# Test Convex Scagnostic
test_that("Convex Scagnotist", {
  expect_gt(sc_convex(data$x,data$y), 0.95)
})
