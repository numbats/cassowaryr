# library(dplyr)
# library(ggplot2)
# library(cassowaryr)
#
# set2 <- anscombe_tidy |>
#   dplyr::filter(set=="2")
#
# calc_scags(set2$x,set2$y)
# plot(set2$x,set2$y)
# set2$x[1:2] <- NA
# calc_scags(set2$x,set2$y)
# plot(set2$x,set2$y)
#
#
# set3 <- anscombe_tidy |>
#   dplyr::filter(set=="3")
#
# calc_scags(set3$x,set3$y)
#
# set4 <- anscombe_tidy |>
#   dplyr::filter(set=="4")
#
# calc_scags(set4$x,set4$y)
#
# ggplot(anscombe_tidy) +
#   geom_point(aes(x=x, y=y))+
#   facet_wrap(~set)
