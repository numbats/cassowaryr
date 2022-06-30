# cassowaryr 2.0.0 (Current)
- In line with the alphahull 2.5 update, cassowaryr is now dependent on interp instead of tripack and has a 
-- `scree()` will throw an error for any scatter plot where there is a large portion of the plot that only contains data in a perfectly straight line (if the convex hull reconnects)
-- `scree()` can take significantly longer on a select few scatter plots
- energy package is now required instead of suggested
- The "line" scatter plot in the `features` data is no longer a perfectly straight line to avoid an error from `scree()` 

# cassowaryr 1.0.1
- All functions now have a value description. 
- The draw functions now have the option for outlier removal. 
- `draw_alphahull()` function fill option now fills the polygon.

# cassowaryr 1.0.0
- This package contains functions to compute scagnostics measuring different patterns in scatter plots.
