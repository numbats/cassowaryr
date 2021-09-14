# Generate some data
library(dplyr)
library(tourr)
library(geozoo)
library(vaast)

set.seed(0903)
df <- tibble(x=1:10, y=1:10, set="line") %>%
  bind_rows(tibble(x=rnorm(100), y=rnorm(100), set="norm"))
d <- as_tibble(sphere.hollow(p=2, n=100)$points) %>%
  rename(x = V1, y = V2) %>%
  mutate(set = "circle")
df <- df %>%
  bind_rows(d)
d <- flea %>%
  select(tars2, head) %>%
  #mutate(head = jitter(head, factor=0.1)) %>%
  rename(x=tars2, y=head) %>%
  mutate(set = "stripes")
df <- df %>%
  bind_rows(d)
d <- flea %>%
  select(tars1, aede1) %>%
  #mutate(tars1 = jitter(tars1, factor=0.1),
  #       aede1 = jitter(aede1, factor=0.1)) %>%
  rename(x=tars1, y=aede1) %>%
  mutate(set = "clumps")
df <- df %>%
  bind_rows(d) %>%
  mutate(set = factor(set,
                      levels = c("line", "norm", "circle", "stripes", "clumps")))

# Compute scagnostics
s <- df %>%
  group_by(set) %>%
  summarise(calc_scags(x, y,
                       scags = c("outlying", "stringy", "striated",
                                 "clumpy", "sparse",
                                 "monotonic", "dcor"))) %>%
  mutate(plot = "") %>%
  select(plot, set, outlying, stringy, striated,
         clumpy, sparse, monotonic, dcor)

