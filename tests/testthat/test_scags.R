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
d1 <- sample(1:10, 100, replace=TRUE)
d2 <- sample(1:10, 100, replace=TRUE)
c1 <- rnorm(100)
c2 <- rnorm(100)

#most striated_adjusted (both discrete)
sc_striated_adjusted(d1,d2)
sc_striated(d1,d2)

#next most (one variable is discrete)
sc_striated_adjusted(d1,c1)
sc_striated(d1,c1)

# least (both continuous)
sc_striated_adjusted(c1,c2)
sc_striated(c1,c2)
