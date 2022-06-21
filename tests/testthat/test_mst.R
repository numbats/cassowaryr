# Load Libraries
library(dplyr)

set.seed(4213)

# Clumpy test

x2 <- c(0,10,10,10,20,80,90,90,90,100,10,90)
y2 <- c(50,60,50,40,50,50,60,50,40,50,100,0)

#plot(x2,y2)

test_that("Clumpy Scagnotist", {
  expect_equal(sc_clumpy(x2,y2), (1-0.4/0.6))
})

# Test for extra clumpy issue
x1 <- c(0,1,4,4.0001, 6, 6.0001, 9, 10)
y1 <- c(10,10,0,0,0,0,10,10)
#sc_clumpy(x1,y1)

# Stringy Test
x1 <- c(0,10,20,30,40,50,60,70,80,90,100, 15)
y1 <- c(0,10,20,30,40,50,60,70,80,90,100, 25)
#plot(x1,y1)

test_that("Stringy Scagnotist", {
  expect_equal(sc_stringy(x1,y1), (8/9))
})

#Striated test
test_that("Striated Scagnotist", {
  expect_equal(sc_striated(x1,y1), (6/9))
})

# Sparse Test
x1 <- c(0,9,17,22,26,26,0,0,0,0,2)
y1 <- c(0,0,0,0,0,2,10,17,23,26,26)
#sample size adjustment
n = 11/500
w = 0.7 + 0.3/(1+n^2)

test_that("Sparse Scagnotist", {
  expect_equal(sc_sparse(x1,y1), w*(9/26), tolerance=0.001)
})

# Skewed Test
test_that("Skewed Scagnotist", {
  expect_equal(sc_skewed(x1,y1), w*(4/7), tolerance=0.001)
})

# Oulying Tests
# No Outliers
x1<- c(0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1)
y1 <- c(0,0,1,2,2,3,4,4,5,6,6,7,8,8,9,10)

# Single outlier
x2<- c(x1, 10)
y2 <- c(y1,10)
#plot(x1,y1)

# Outlier with a Close Value
x3 <- c(x2, 10)
y3 <- c(y2, 9)

# two connected outliers
x4 <- c(x2, 10)
y4 <- c(y2, 7)

# internal outlier
x5 <- c(x4, 10)
y5 <- c(y4, 6)

# internal outlier
x6 <- c(0,1,1,1,0,0,0,0,0,0,0,1,1,1, 10, 10, 10)
y6 <- c(0,0,1,2,2,3,4,5,6,7,8,8,9,10, 10, 6, 7)

#plot(x2,y2)

test_that("Outlying Scagnotist", {
  expect_equal(sc_outlying(x1,y1), 0) #no outliers
  expect_equal(sc_outlying(x2,y2), (9/24)) #single external outleirs
  expect_equal(sc_outlying(x3,y3), 0) #close edge to outlier so no longer outlier
  expect_equal(sc_outlying(x4,y4), 12/27) # internal and external outliers connected
  expect_equal(sc_outlying(x5,y5), 3/28) #single internal outlier 1 (becomes single external outlier)
  expect_equal(sc_outlying(x6,y6), 12/26) #single internal outlier 2
})

# Striated2 Test
x1 <- c(0,5,0,10)
y1 <- c(0,5,10,10)

#slightly more complicated
x2 <- c(0,5,0,10,10, 2.5, 0)
y2 <- c(0,5,10,10,0, 5, 5)

test_that("Striated2 Tests", {
  expect_equal(sc_striated2(x1,y1), 1)
  expect_equal(sc_striated2(x2,y2), 0.8)
})

# Triangle collapsed data
x <- seq(15)
y <- c(8,9,10,9,8,7,6,5,4,3,2,1,0,1,2)
