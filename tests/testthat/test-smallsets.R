# PLACEHOLDER FOR CHECKS
# Remove missing values
d <- tibble::tibble(x=x, y=y)
d <- d[stats::complete.cases(d),]
if (length(x) > nrow(d))
  message("WARNING: ", length(x)-nrow(d), " observations in have been removed. \n")
x <- d$x
y <- d$y

if(is.null(sm_list)){
  # this is null when one of the variables is constant after outlier removal
  return(
    dplyr::tibble("outlying"=NA,
                  "stringy"=NA,
                  "striated"=NA,
                  "grid" = NA,
                  "clumpy"=NA,
                  "clumpy2"=NA,
                  "sparse"=NA,
                  "skewed"=NA,
                  "convex"=NA,
                  "skinny"=NA,
                  "monotonic"=NA,
                  "splines"=NA,
                  "dcor"=NA,
                  "stripes"=NA)
  )
}

# require(ggplot2)
# require(tidyr)
# require(dplyr)
#
# data(anscombe_tidy)
# ggplot(anscombe_tidy, aes(x=x, y=y)) +
#   geom_point() +
#   facet_wrap(~set, ncol=2, scales = "free")
#
# draw_mst(anscombe$x1, anscombe$y1)
# draw_mst(anscombe$x2, anscombe$y2, out.rm=TRUE)
# draw_mst(anscombe$x3, anscombe$y3, out.rm=TRUE)
# draw_mst(anscombe$x4, anscombe$y4, out.rm=TRUE)

# sc_stringy(anscombe$x1, anscombe$y1)
# sc_stringy(anscombe$x2, anscombe$y2)
# sc_stringy(anscombe$x3, anscombe$y3)
# sc_stringy(anscombe$x4, anscombe$y4)
#
# sc_striated(anscombe$x1, anscombe$y1)
# sc_striated(anscombe$x2, anscombe$y2)
# sc_striated(anscombe$x3, anscombe$y3)
# sc_striated(anscombe$x4, anscombe$y4)

# test_that("multiplication works", {
#   expect_equal(2 * 2, 4)
# })
