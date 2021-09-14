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

