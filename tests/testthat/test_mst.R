# Load Libraries
library(dplyr)

set.seed(4213)

# Clumpy test

# Case 1
x1 <- c(0,10,10,10,20,80,90,90,90,100,10,90)
y1 <- c(50,60,50,40,50,50,60,50,40,50,100,0)

# Case 2 (similar but without long edge)
x2 <- c(0,1,1,1,2, 8,9,9,9,10)
y2 <- c(1,2,1,0,1, 9,10,9,8,9)

e_length <- sqrt((0.7)^2 + (0.7)^2)

test_that("Clumpy Scagnostic", {
  expect_equal(sc_clumpy(x1,y1), (1-0.4/0.6))
  expect_equal(sc_clumpy(x2,y2), (1-0.1/e_length))
})



# Test for extra clumpy issue
# x1 <- c(0,1,4,4.0001, 6, 6.0001, 9, 10)
# y1 <- c(10,10,0,0,0,0,10,10)
# sc_clumpy(x1,y1, binner=NULL)
# sc_clumpy(x1,y1)

# Stringy Test
x1 <- c(0,10,20,30,40,50,60,70,80,90,100, 15)
y1 <- c(0,10,20,30,40,50,60,70,80,90,100, 25)
#plot(x1,y1)

test_that("Stringy Scagnostic", {
  expect_equal(sc_stringy(x1,y1), (8/9))
  expect_equal(sc_stringy2(x1,y1), 1)
})


# additional tests for stringy index
test_that("Stringy scagnostic", {
  set.seed(1200)

  x  <- seq(0, 1, length.out = 50)

  y1 <- x
  y2 <- 2 * x

  s1 <- suppressWarnings(sc_stringy(x, y1))
  s2 <- suppressWarnings(sc_stringy(x, y2))
  s3 <- suppressWarnings(sc_stringy2(x, y1))
  s4 <- suppressWarnings(sc_stringy2(x, y1))

  # we expect exactly 1
  expect_equal(s1, 1, tolerance = 1e-8)
  expect_equal(s2, 1, tolerance = 1e-8)
  expect_equal(s3, 1, tolerance = 1e-8)
  expect_equal(s4, 1, tolerance = 1e-8)
})

#Striated test
test_that("Striated Scagnostic", {
  expect_equal(sc_striated(x1,y1), (6/9))
})

# Sparse Test
x1 <- c(0,9,17,22,26,26,0,0,0,0,2)
y1 <- c(0,0,0,0,0,2,10,17,23,26,26)
#sample size adjustment
n = 11/500
w = 0.7 + 0.3/(1+n^2)

test_that("Sparse Scagnostic", {
  expect_equal(sc_sparse(x1,y1), w*(9/26), tolerance=0.001)
})

# Skewed Test
test_that("Skewed Scagnostic", {
  expect_equal(sc_skewed(x1,y1), w*(4/7), tolerance=0.001)
})

# Oulying Tests
# No Outliers
x1<- c(0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1)
y1 <- c(0,0,1,2,2,3,4,4,5,6,6,7,8,8,9,10)

# Two outliers
x2 <- c(x1, 10)
y2 <- c(y1,10)
# plot(x2,y2)
# draw_(x2,y2)

# Outlier with a Close Value
x3 <- c(x2, 10)
y3 <- c(y2, 9)
# plot(x3,y3)
# draw_mst(x3,y3)

# two connected outliers
x4 <- c(x2, 10)
y4 <- c(y2, 7)
# plot(x4,y4)
# draw_mst(x4,y4)

# internal outlier
x5 <- c(x4, 10)
y5 <- c(y4, 6)
# plot(x5,y5)
# draw_mst(x5,y5)

# internal outlier
x6 <- c(0,1,1,1,0,0,0,0,0,0,0,1,1,1, 10, 10, 10)
y6 <- c(0,0,1,2,2,3,4,5,6,7,8,8,9,10, 10, 6, 7)
# plot(x6,y6)
# draw_mst(x6,y6)

#plot(x2,y2)

test_that("Outlying Scagnostic", {
  expect_equal(sc_outlying(x1,y1), 0) #no outliers
  expect_equal(sc_outlying(x2,y2), (9/24)) #single external outleirs
  expect_equal(sc_outlying(x3,y3), 0) #close edge to outlier so no longer outlier
  expect_equal(sc_outlying(x4,y4), 12/27) # internal and external outliers connected
  expect_equal(sc_outlying(x5,y5), 3/28) #single internal outlier 1 (becomes single external outlier)
  expect_equal(sc_outlying(x6,y6), 12/26) #single internal outlier 2
})

# grid Test
x1 <- c(0,5,0,10)
y1 <- c(0,5,10,10)

#slightly more complicated
x2 <- c(0,5,0,10,10, 2.5, 0)
y2 <- c(0,5,10,10,0, 5, 5)

test_that("grid Tests", {
  expect_equal(sc_grid(x1,y1), 1)
  expect_equal(sc_grid(x2,y2), 0.8)
  expect_warning(calc_scags(x2,y2, "striated2"))
})

# Triangle collapsed data
x <- seq(15)
y <- c(8,9,10,9,8,7,6,5,4,3,2,1,0,1,2)
