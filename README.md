
<!-- README.md is generated from README.Rmd. Please edit that file -->

# About

<!-- badges: start -->
<!-- badges: end -->

The `vaast` package provides functions to compute scagnostics on pairs
of numeric variables in a data set.

The term **scagnostics** refers to scatter plot diagnostics. This is a
collection of techniques for automatically extracting interesting visual
features from pairs of variables. This package is an implementation of
graph theoretic scagnostics developed by Wilkinson, Anand, and Grossman
(2005) in pure R.

## Installation

``` r
remotes::install_github("numbats/vaast")
```

## Examples

``` r
data("anscombe_tidy")
sc_outlying(anscombe$x1, anscombe$y1)
sc_outlying(anscombe$x2, anscombe$y2)
sc_outlying(anscombe$x3, anscombe$y3)
sc_outlying(anscombe$x4, anscombe$y4)
```

## About the calculations

### Graph-based measures

A 2-d scatter plot can be represented by a combination of three graphs
which are computed directly from the Delauney-Voroni tesselation.

1.  A **minimum spanning tree** weighted by the lengths of the Delauney
    triangles
2.  The **convex hull** of the points i.e. the outer segments of the
    triangulation
3.  The **alpha hull ** (also called concave hull) i.e. formed by
    connect the outer edges of triangles that are enclosed within a ball
    of radius *alpha*.

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
