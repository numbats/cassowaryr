# Basic test case
x <- c(0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1)
y <- c(0,0,1,2,2,3,4,4,5,6,6,7,8,8,9,10)

# Data pre-processing
xy <- cbind(unitize(x), unitize(y))
del <- alphahull::delvor(xy)
weights <- gen_edge_lengths(del)
mst <- gen_mst(del, weights)
mst_weights <- igraph::E(mst)$weight

# Manual value calculations for each alpha
manual_weights <- c(rep(0.1, 10), rep(0.2, 4), 1)
manual_q90 <- 0.2
manual_rahman <- sqrt(sum(c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.1))/16)
manual_omega <- 0.2 + 1.5 * 0.1

# run tests
test_that("Test alpha calculations", {
  # q90
  expect_equal(alpha_q90(mst_weights), manual_q90)

  # rahman
  expect_equal(alpha_rahman(mst_weights), manual_rahman)

  # omega
  expect_equal(alpha_omega(mst_weights), manual_omega)
})


test_that("Test get_numeric_alpha", {
  # q90
  expect_equal(alpha_q90(mst_weights),
               get_numeric_alpha("q90", del, weights))

  # rahman
  expect_equal(alpha_rahman(mst_weights),
               get_numeric_alpha("rahman", del, weights))

  # omega
  expect_equal(alpha_omega(mst_weights),
               get_numeric_alpha("omega", del, weights))
})


