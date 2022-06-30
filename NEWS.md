# cassowaryr 2.0.0 (Current)

- In line with the alphahull 2.5 update, cassowaryr is now dependent on the package interp instead of tripack and which results in some issues to be addressed: 
-- `scree()` which calculates the convex and alpha hulls will return an error for any scatter plot where most of the data lies on a straight line 
-- `scree()` can take significantly longer on a select few scatter plots
- The energy package is now required instead of suggested
- The "line" subset in the `features` data is no longer a perfectly straight line to avoid the error above.

# cassowaryr 1.0.1

- All functions now have a value description. 
- The draw functions now have the option for outlier removal. 
- `draw_alphahull()` function fill option now fills the polygon.

# cassowaryr 1.0.0

- This package contains functions to compute scagnostics measuring different patterns in scatter plots.
