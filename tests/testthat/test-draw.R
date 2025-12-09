library(cassowaryr)
library(ggplot2)

# No Outliers
x1 <- c(0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1)
y1 <- c(0,0,1,2,2,3,4,4,5,6,6,7,8,8,9,10)
# plot of points
ggplot() + geom_point(aes(x1,y1))
# no outlier removal
draw_mst(x1,y1, outlier_rm = FALSE) +
  xlim(0,1) + ylim(0,1)
# outlier removal, old version
draw_mst_old(x1,y1, out.rm=TRUE)+
  xlim(0,1) + ylim(0,1)
# new version
draw_mst(x1,y1, outlier_rm = TRUE)+
  xlim(0,1) + ylim(0,1)


# Two outliers
x2 <- c(x1, 10)
y2 <- c(y1,10)
# plot of points
ggplot() + geom_point(aes(x2,y2))
# no outlier removal
draw_mst(x2,y2, outlier_rm = FALSE) +
  xlim(0,1) + ylim(0,1)
# outlier removal, old version
draw_mst_old(x2,y2, out.rm=TRUE)+
  xlim(0,1) + ylim(0,1)
# new version
draw_mst(x2,y2, outlier_rm = TRUE)+
  xlim(0,1) + ylim(0,1)


# Outlier with a Close Value
x3 <- c(x2, 10)
y3 <- c(y2, 9)
# plot of points
ggplot() + geom_point(aes(x3,y3))
# no outlier removal
draw_mst(x3,y3, outlier_rm = FALSE) +
  xlim(0,1) + ylim(0,1)
# outlier removal, old version
draw_mst_old(x3,y3, out.rm=TRUE)+
  xlim(0,1) + ylim(0,1)
# new version
draw_mst(x3,y3, outlier_rm = TRUE)+
  xlim(0,1) + ylim(0,1)

# two connected outliers
x4 <- c(x2, 10)
y4 <- c(y2, 7)
# plot of points
ggplot() + geom_point(aes(x4,y4))
# no outlier removal
draw_mst(x4,y4, outlier_rm = FALSE) +
  xlim(0,1) + ylim(0,1)
# outlier removal, old version
draw_mst_old(x4,y4, out.rm=TRUE)+
  xlim(0,1) + ylim(0,1)
# new version
draw_mst(x4,y4, outlier_rm = TRUE)+
  xlim(0,1) + ylim(0,1)

# internal outlier
x5 <- c(x4, 10)
y5 <- c(y4, 6)
# plot of points
ggplot() + geom_point(aes(x5,y5))
# no outlier removal
draw_mst(x5,y5, outlier_rm = FALSE) +
  xlim(0,1) + ylim(0,1)
# outlier removal, old version
draw_mst_old(x5,y5, out.rm=TRUE)+
  xlim(0,1) + ylim(0,1)
# new version
draw_mst(x5,y5, outlier_rm = TRUE)+
  xlim(0,1) + ylim(0,1)

# internal outlier
x6 <- c(0,1,1,1,0,0,0,0,0,0,0,1,1,1, 10, 10, 10)
y6 <- c(0,0,1,2,2,3,4,5,6,7,8,8,9,10, 10, 6, 7)
# plot of points
ggplot() + geom_point(aes(x6,y6))
# no outlier removal
draw_mst(x6,y6, outlier_rm = FALSE) +
  xlim(0,1) + ylim(0,1)
# outlier removal, old version
draw_mst_old(x6,y6, out.rm=TRUE)+
  xlim(0,1) + ylim(0,1)
# new version
draw_mst(x6,y6, outlier_rm = TRUE)+
  xlim(0,1) + ylim(0,1)
