library(cassowaryr)
library(dplyr)
library(alphahull)
library(ggplot2)
library(vdiffr)

cl <- features |> filter(feature == "clusters")

test_that("draw function tests", {

  p1 <- draw_alphahull(cl$x, cl$y)
  expect_doppelganger("draw_alphahull", p1)

  p2 <- draw_mst(cl$x, cl$y)
  expect_doppelganger("draw_mst", p2)

  p3 <- draw_convexhull(cl$x, cl$y)
  expect_doppelganger("draw_convexhull", p3)
}
)
