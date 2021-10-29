#generate sin wave test
x <- seq(-pi, pi, 0.01)
y <- sin(x)

test_that("dcor vs splines", {
  expect_false(isTRUE(all.equal(sc_splines(x,y), sc_dcor(x,y))))
})
