
<!-- README.md is generated from README.Rmd. Please edit that file -->

# About

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/numbats/cassowaryr/branch/master/graph/badge.svg)](https://app.codecov.io/gh/numbats/cassowaryr?branch=master)
[![R-CMD-check](https://github.com/numbats/cassowaryr/workflows/R-CMD-check/badge.svg)](https://github.com/numbats/cassowaryr/actions)
<!-- badges: end -->

The `cassowaryr` package provides functions to compute scagnostics on
pairs of numeric variables in a data set.

The term **scagnostics** refers to scatter plot diagnostics, originally
described by John and Paul Tukey. This is a collection of techniques for
automatically extracting interesting visual features from pairs of
variables. This package is an implementation of graph theoretic
scagnostics developed by Wilkinson, Anand, and Grossman (2005) in pure
R.

## Installation

The package can be installed from CRAN using

> `install.packages("cassowaryr")`

and from GitHub using

> `remotes::install_github("numbats/cassowaryr")`

to install the development version.

## Examples

``` r
library(cassowaryr)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
# A single scagnostic on two vectors
data("anscombe_tidy")
sc_outlying(anscombe$x1, anscombe$y1)
#> [1] 0
```

``` r
data("datasaurus_dozen")
datasaurus_dozen %>%
  dplyr::group_by(dataset)%>%
  dplyr::summarise(calc_scags(x, y, scags=c("clumpy2", "monotonic")))
#> Registered S3 method overwritten by 'cli':
#>   method     from         
#>   print.boxx spatstat.geom
#> # A tibble: 13 x 3
#>    dataset    clumpy2 monotonic
#>    <chr>        <dbl>     <dbl>
#>  1 away        0        0.0573 
#>  2 bullseye    0.811    0.0787 
#>  3 circle      0        0.0773 
#>  4 dino        0.0122   0.0651 
#>  5 dots        0.993    0.126  
#>  6 h_lines     0.912    0.0520 
#>  7 high_lines  0.914    0.00287
#>  8 slant_down  0.888    0.0669 
#>  9 slant_up    0.862    0.0861 
#> 10 star        0.919    0.0514 
#> 11 v_lines     0.941    0.0566 
#> 12 wide_lines  0.920    0.0522 
#> 13 x_shape     0.882    0.0205
```

## About the name

CAlculate Scagnostics on Scatterplots Over Wads of Associated Real
numberYs in R

## About the calculations

### Graph-based measures

A 2-d scatter plot can be represented by a combination of three graphs
which are computed directly from the Delauney-Voroni tesselation.

1.  A **minimum spanning tree** weighted by the lengths of the Delauney
    triangles
2.  The **convex hull** of the points i.e. the outer segments of the
    triangulation
3.  The **alpha hull** (also called concave hull) i.e. formed by connect
    the outer edges of triangles that are enclosed within a ball of
    radius *alpha*.

All graph based scagnostic measures are computed with respect to these
three graphs.

Prior to graph construction decisions must be made about filtering
outliers (done with respect to the distribution of edge lengths in the
triangulation) and thinning the size of the graphs by performing binning
(for computational speed). For the moment we can forge ahead without
these but it is worth keeping in mind that the package needs to be
flexible enough to include them.  
There are also opportunities to experiment with the preprocessing here.

Two MST measures “clumpy” and “outlying” are known to cause problems.

As all the graph based measures rely on the triangulation, they could be
computed lazily. More concretely, if you are only interested in
computing “skinny” you don’t need to compute the spanning tree. If you
computed “skinny” but wanted to then compute “convex” you shouldn’t need
to reconstruct the alpha-hull and so on… To begin let’s not worry about
that and focus on implementations of each measure.

### Association-based measures

These are computed directly from the 2-d point clouds, and do not need
to be constructed from the graph.
