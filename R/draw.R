# Code for drawing alphahull,
# not in package format yet

library(geozoo)
library(ggplot2)
library(vaast)
library(alphahull)
set.seed(0903)
d <- as_tibble(sphere.hollow(p=2, n=300)$points) %>%
rename(x = V1, y = V2) %>%
filter(!(x > -0.5 & y > -0.5)) %>%
mutate(x = x + rnorm(length(x), sd=0.2), y = y + rnorm(length(y), sd=0.2))
p <- ggplot(d, aes(x, y)) + geom_point()

d_ahull <- ahull(d$x, d$y, a=0.5)
d_ahull_c <- d_ahull$ashape.obj
p + geom_segment(data=as_tibble(d_ahull_c$edges), aes(x=x1, xend=x2, y=y1, yend=y2))

