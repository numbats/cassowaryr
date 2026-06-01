# Unresolved issues: 1, 6, 7

# Issue 1: "triangle collapse" & "points are colinear" warning from clumpy test
x2 <- c(0,1,1,1,2, 8,9,9,9,10)
y2 <- c(1,2,1,0,1, 9,10,9,8,9)
sc_clumpy(x2,y2)

# Cause: Data itself, not the scagnostic.
# The "rahman" alphahull value is making janky alphahull, see
draw_alphahull(x2,y2) # warning
draw_alphahull(x2,y2, alpha = "omega") # fine

# would fix with
# sc_clumpy(x2,y2, alpha = "omega")
# but clumpy cannot control alpha (obviously because it is not used)
# confusing that the alphahull is throwing a warning for mst measure

# Open an issue on this in the github.
# For now, I should just change the data that goes in the test.

# Issue 2: Stringy scagnostic values different to expected
x1 <- c(0,10,20,30,40,50,60,70,80,90,100, 15)
y1 <- c(0,10,20,30,40,50,60,70,80,90,100, 25)
sc_stringy05(x1,y1)

# Cause: I have no idea why the value was set to be 1.
# that would have never been correct for this data?
# Should be: longest path/total mst length
path <- 10*10*sqrt(2)
extra <- sqrt(sum((c(20,20) - c(15,25))^2))
total <- path +  extra
path/total

# Issue 3: Grid scagnostic giving different values than expected
x2 <- c(0,5,0,10,10, 2.5, 0)
y2 <- c(0,5,10,10,0, 5, 5)
sc_grid(x2,y2)

# The solutions say this should be 0.8, but I have no idea why. 0.6 looks like the correct answer to me.

# Issue 4: Scree & binningstests not working
# Was calling wrong file path

# Issue 5: Scree tests assumed old defaults
# Fixed by swapping defaults


# Issue 6: Skinny test
# We get perfect colinear warning from
sc_skinny(x, y1) # since straight line
# but not
sc_skinny(x, y2) # since it is a curved line
# I should get rid of the perfectly colinear warning
# & also ask Tina about it. She had this set to 0.7, but it gives a 1.


# Issue 7 : Tina Skinny Code
# Ask Tina what this code is supposed to do.
# It blows the computation time of skinny out so much it can't really be used
# If it is checking for something in the alphahull, it should be checked when
# we compute the alphahull not

# sc <- scree(x, y, alpha = alpha, out.rm = out.rm, binner = binner)
# if (is.null(sc$del)) {
#   return(1)
# }
#
# dc <- sc_dcor(x, y)
# d2 <- mst_deg2(x, y)
#
# if (dc >= 1 - 1e-8 || (!is.na(d2) && d2 > 0.9)) {
#   return(skinny_fix(x, y))
# }
#
# sc_skinny.scree(sc)
