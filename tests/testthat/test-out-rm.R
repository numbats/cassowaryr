
# current test case
x <- c(x1, 10)
y <- c(y1,10)

xy <- cbind(unitize(x), unitize(y))
del <- alphahull::delvor(xy)
weights <- gen_edge_lengths(del)
mst <- gen_mst(del, weights)
mst_weights <- igraph::E(mst)$weight

# TEST CASES
# No Outliers
x1 <- c(0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1)
y1 <- c(0,0,1,2,2,3,4,4,5,6,6,7,8,8,9,10)

xy1 <- cbind(unitize(x1), unitize(y1))
del1 <- alphahull::delvor(xy1)
weights1 <- gen_edge_lengths(del1)
mst1 <- gen_mst(del1, weights1)
mst_weights1 <- igraph::E(mst1)$weight

# One outlier
x2 <- c(x1, 10)
y2 <- c(y1,10)

xy2 <- cbind(unitize(x2), unitize(y2))
del2 <- alphahull::delvor(xy2)
weights2 <- gen_edge_lengths(del2)
mst2 <- gen_mst(del2, weights2)
mst_weights2 <- igraph::E(mst2)$weight

# Outlier with a Close Value (so no outliers)
x3 <- c(x2, 10)
y3 <- c(y2, 9)

xy3 <- cbind(unitize(x3), unitize(y3))
del3 <- alphahull::delvor(xy3)
weights3 <- gen_edge_lengths(del3)
mst3 <- gen_mst(del3, weights3)
mst_weights3 <- igraph::E(mst3)$weight

# Two connected outliers
x4 <- c(x2, 10)
y4 <- c(y2, 7)

xy4 <- cbind(unitize(x4), unitize(y4))
del4 <- alphahull::delvor(xy4)
weights4 <- gen_edge_lengths(del4)
mst4 <- gen_mst(del4, weights4)
mst_weights4 <- igraph::E(mst4)$weight


# Internal outlier
x5 <- c(0,1,1,1,0,0,0,0,0,0,0,1,1,1, 10, 10, 10)
y5 <- c(0,0,1,2,2,3,4,5,6,7,8,8,9,10, 10, 6, 7)

xy5 <- cbind(unitize(x5), unitize(y5))
del5 <- alphahull::delvor(xy5)
weights5 <- gen_edge_lengths(del5)
mst5 <- gen_mst(del5, weights5)
mst_weights5 <- igraph::E(mst5)$weight


# Iterative outlier
x6 <- c(0,1,1,2,2,3,3,4,4,5,5,6,7.2,8.5,10)
y6 <- c(0,0,1,1,2,2,3,3,4,4,5,6,7.2,8.5,10)

xy6 <- cbind(unitize(x6), unitize(y6))
del6 <- alphahull::delvor(xy6)
weights6 <- gen_edge_lengths(del6)
mst6 <- gen_mst(del6, weights6)
mst_weights6 <- igraph::E(mst6)$weight

test_that("Test outlying_identify", {
  # no outliers
  expect_equal(outlying_identify(mst1, mst_weights1), NULL)
  # one edge outlier
  expect_equal(outlying_identify(mst2, mst_weights2), 17)
  # long edge next to close edge (not outlier)
  expect_equal(outlying_identify(mst3, mst_weights3), NULL)
  # two connected outliers
  expect_equal(outlying_identify(mst4, mst_weights4), c(17, 18))
  # internal outlier
  expect_equal(outlying_identify(mst5, mst_weights5), 15)
  # iterative outlier (only gets first round)
  expect_equal(outlying_identify(mst6, mst_weights6), c(14, 15))
})

test_that("Test outlier_removal", {
  # no outliers
  expect_equal(outlier_removal(xy1), alphahull::delvor(xy1))
  # one edge outlier
  expect_equal(outlier_removal(xy2), alphahull::delvor(xy2[-17,]))
  # long edge next to close edge (not outlier)
  expect_equal(outlier_removal(xy3), alphahull::delvor(xy3))
  # two connected outliers
  expect_equal(outlier_removal(xy4), alphahull::delvor(xy4[-c(17, 18),]))
  # internal outlier
  expect_equal(outlier_removal(xy5), alphahull::delvor(xy5[-15,]))
  # iterative outliers
  expect_equal(outlier_removal(xy6), alphahull::delvor(xy6[-c(12:15),]))
})

