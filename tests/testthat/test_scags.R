# Generate some data
library(tidyr)
library(dplyr)
library(tourr)
library(geozoo)
library(vaast)
library(GGally)

# Alternative data set
set.seed(55555)
d_trend <- tibble(x = runif(100) - 0.5) %>%
  mutate(positive=4*x+rnorm(100)*0.5,
         none=rnorm(100)*0.5,
         negative=-4*x+rnorm(100)*0.5) %>%
  pivot_longer(cols=positive:negative, names_to="trend", values_to="y") %>%
  mutate(feature = factor(trend,
                        levels=c("positive", "none", "negative"))) %>%
  select(feature, x, y) %>%
  filter(feature == "positive")

d_strength <- tibble(x = runif(100) -0.5) %>%
  mutate(strong=4*x+rnorm(100)*0.5,
         moderate=4*x+rnorm(100),
         weak=-4*x+rnorm(100)*3) %>%
  pivot_longer(cols=strong:weak,
               names_to="strength", values_to="y") %>%
  mutate(feature = factor(strength,
                           levels=c("strong", "moderate", "weak"))) %>%
  select(feature, x, y) %>%
  filter(feature == "weak")

d_form <- tibble(x = runif(100) -0.5) %>%
  mutate(linear=4*x+rnorm(100)*0.5,
         nonlinear1=12*x^2+rnorm(100)*0.5,
         nonlinear2=2*x - 5*x^2 +rnorm(100)*0.1) %>%
  pivot_longer(cols=linear:nonlinear2, names_to="feature",
               values_to="y") %>%
  select(feature, x, y) %>%
  filter(feature == "nonlinear2")

d_outliers <- tibble(x = runif(100) -0.5) %>%
  mutate(y=4*x+rnorm(100)*0.5)
d_outliers <- d_outliers %>%
  bind_rows(tibble(x=runif(5)/10-0.45, y=2+rnorm(5)*0.5)) %>%
  mutate(feature = "outliers") %>%
  select(feature, x, y)

d_clusters <- tibble(x = c(rnorm(50)/6 - 0.5,
                           rnorm(50)/6,
                           rnorm(50)/6 + 0.5)) %>%
  mutate(y=c(rnorm(50)/6,
             rnorm(50)/6+1, rnorm(50)/6)) %>%
  mutate(feature = "clusters") %>%
  select(feature, x, y)

d_gaps <- tibble(x = runif(150)) %>%
  mutate(y=runif(150))
d_gaps <- d_gaps %>%
  filter(!(between(x+2*y, 1.2, 1.6))) %>%
  mutate(feature = "gaps") %>%
  select(feature, x, y)

d_barrier <- tibble(x = runif(200)) %>%
  mutate(y=runif(200))
d_barrier <- d_barrier %>%
  filter(-x+3*y<1.2) %>%
  mutate(feature = "barrier") %>%
  select(feature, x, y)

l_shape <- tibble(x = c(rexp(50, 0.01), runif(50)*20),
                  y = c(runif(50)*20, rexp(50, 0.01))) %>%
  mutate(feature = "l-shape") %>%
  select(feature, x, y)

discrete <- tibble(x = rnorm(200)) %>%
  mutate(y = -x+rnorm(25)*0.1 + rep(0:7, 25)) %>%
  filter((scale(x)^2 + scale(y)^2) < 2) %>%
  mutate(feature = "discrete") %>%
  select(feature, x, y)

hetero <- tibble(x = runif(200)-0.5) %>%
  mutate(y = -2*x+rnorm(200)*(x+0.5)) %>%
  mutate(feature = "heteroskedastic") %>%
  select(feature, x, y)

features <- bind_rows(d_trend,
                      d_strength,
                      d_form,
                      d_outliers,
                      d_clusters,
                      d_gaps,
                      d_barrier,
                      discrete,
                      l_shape)
s <- features %>%
  group_by(feature) %>%
  summarise(calc_scags(x, y,
                       scags = c("outlying", "stringy", "striated",
                                 "clumpy", "sparse",
                                 "monotonic", "dcor")))
ggpairs(s, columns = 2:8)

