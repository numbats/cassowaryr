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
