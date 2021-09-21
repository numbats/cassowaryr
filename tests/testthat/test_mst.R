# Load Libraries
library(dplyr)
library(vaast)
library(usethis)

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
sc_clumpy(x1,y1)

# Stringy Test
x1 <- c(0,10,20,30,40,50,60,70,80,90,100, 15)
y1 <- c(0,10,20,30,40,50,60,70,80,90,100, 25)
#plot(x1,y1)

test_that("Stringy Scagnotist", {
  expect_equal(sc_stringy(x1,y1), (8/9))
})

#Striated test
x1 <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)
y1 <- c(8,9,10,9,8,7,6,5,4,3,2,1,0,1,2)
#plot(x1,y1)

test_that("Striated Scagnotist", {
  expect_equal(sc_striated(x1,y1), (11/15))
})

# Sparse Test
x1 <- c(0,9,17,22,26,26,0,0,0,0,2)
y1 <- c(0,0,0,0,0,2,10,17,23,26,26)
#plot(x1,y1)
test_that("Sparse Scagnotist", {
  expect_equal(sc_sparse(x1,y1), (9/26))
})

# Skewed Test
#plot(x1,y1)
test_that("Sparse Scagnotist", {
  expect_equal(sc_skewed(x1,y1), (4/7))
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
x3 <- c(x1, 10)
y3 <- c(y1, 9)
#plot(x2,y2)

test_that("Outlying Scagnotist", {
  expect_equal(sc_outlying(x1,y1), 0)
  expect_equal(sc_outlying(x2,y2), (9/24))
  expect_equal(sc_outlying(x3,y3), (10/25))
})
