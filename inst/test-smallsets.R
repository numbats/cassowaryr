# PLACEHOLDER FOR CHECKS

# ISSUE 1
# If the data is a perfectly straight line BEFORE outlier removal, we get an error
x  <- seq(0, 1, length.out = 50)
draw_mst(x, x)

# If the data is a perfectly straight line AFTER outlier removal, the MST will work
draw_mst(anscombe$x3, anscombe$y3, out.rm=TRUE)

# HOWEVER if I remove the check for straight lines that comes before outlier removal
# and then rerun
draw_mst(x, x)
# I get "Error: shull: linear structure, aborting"
# Which is the same error I get with a perfectly horizontal line

# OK, it appears that the anscombe$x3 data is NOT actually a PERFECTLY straight line.
# This is the source of the error.

# Solutions: 1) add a tiny bit of error (already happens if triangle collapses, so might as well?)

# Check all work
sc_skinny(x, x)
sc_stringy06(x, x)

# NEW ISSUE
draw_mst(anscombe$x4, anscombe$y4, out.rm=TRUE)
# when plotting the MST for anscomb, it doesnt change.
# should set maximum and minimum as 0 and 1 before plotting?
# Gives better idea of structure anyways


# NEW ISSUE
# ??? These striated values make no sense
sc_striated(anscombe$x2, anscombe$y2)
sc_striated(anscombe$x3, anscombe$y3)
sc_striated(anscombe$x4, anscombe$y4)

# look at the MST
draw_mst(anscombe$x1, anscombe$y1)
draw_mst(anscombe$x2, anscombe$y2, out.rm=TRUE)
draw_mst(anscombe$x3, anscombe$y3, out.rm=TRUE)


require(ggplot2)
require(tidyr)
require(dplyr)

data(anscombe_tidy)
ggplot(anscombe_tidy, aes(x=x, y=y)) +
  geom_point() +
  facet_wrap(~set, ncol=2, scales = "free")

#
draw_mst(anscombe$x1, anscombe$y1)
draw_mst(anscombe$x2, anscombe$y2, out.rm=TRUE)
draw_mst(anscombe$x3, anscombe$y3, out.rm=TRUE)
draw_mst(anscombe$x4, anscombe$y4, out.rm=TRUE)

draw_alphahull(anscombe$x1, anscombe$y1)
draw_alphahull(anscombe$x2, anscombe$y2, out.rm=TRUE)
draw_alphahull(anscombe$x3, anscombe$y3, out.rm=TRUE)
draw_alphahull(anscombe$x4, anscombe$y4, out.rm=TRUE)

sc_skinny(anscombe$x1, anscombe$y1)
sc_skinny(anscombe$x2, anscombe$y2)
sc_skinny(anscombe$x3, anscombe$y3)
sc_skinny(anscombe$x4, anscombe$y4)

sc_stringy06(anscombe$x1, anscombe$y1)
sc_stringy06(anscombe$x2, anscombe$y2)
sc_stringy06(anscombe$x3, anscombe$y3)
sc_stringy06(anscombe$x4, anscombe$y4)
#
sc_striated(anscombe$x1, anscombe$y1)
sc_striated(anscombe$x2, anscombe$y2)
sc_striated(anscombe$x3, anscombe$y3)
sc_striated(anscombe$x4, anscombe$y4)

# lines work when



# IDK what this is
# Remove missing values
d <- tibble::tibble(x=x, y=y)
d <- d[stats::complete.cases(d),]
if (length(x) > nrow(d))
  message("WARNING: ", length(x)-nrow(d), " observations in have been removed. \n")
x <- d$x
y <- d$y

if(is.null(sm_list)){
  # this is null when one of the variables is constant after outlier removal
  return(
    dplyr::tibble("outlying"=NA,
                  "stringy"=NA,
                  "striated"=NA,
                  "grid" = NA,
                  "clumpy"=NA,
                  "clumpy2"=NA,
                  "sparse"=NA,
                  "skewed"=NA,
                  "convex"=NA,
                  "skinny"=NA,
                  "monotonic"=NA,
                  "splines"=NA,
                  "dcor"=NA,
                  "stripes"=NA)
  )
}

