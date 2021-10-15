# Generate some data
library(dplyr)
library(cassowaryr)
library(GGally)

s <- features %>%
  group_by(feature) %>%
  summarise(calc_scags(x, y,
                       scags = c("outlying", "stringy",
                                 "striated",
                                 "clumpy", "sparse",
                                 "monotonic", "dcor")))
ggpairs(s, columns = 2:8)


# Compare striated to striated adjusted

# Simulated data for what it is looking for
d1 <- sample(1:10, 1000, replace=TRUE)
d2 <- sample(1:10, 1000, replace=TRUE)
c1 <- rnorm(1000)
c2 <- rnorm(1000)
#rotated
x3 <- c(c1, c2)
y3 <- c(1.5*c1+1, 1.5*c2+3)

#most striated_adjusted (both discrete)
sc_striated2(d1,d2)
sc_striated(d1,d2)
plot(d1,d2)

#next most (one variable is discrete)
sc_striated2(d1,c1)
sc_striated(d1,c1)
plot(d1,c1)

# least (both continuous)
sc_striated2(c1,c2)
sc_striated(c1,c2)
plot(c1,c2)

# rotated (both continuous)
sc_striated2(x3,y3)
sc_striated(x3,y3)
plot(x3,y3)

