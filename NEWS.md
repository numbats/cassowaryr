# cassowaryr 2.0.11
- Broke scree function up so binning, outlier removal, and alpha calculations are in their own functions (and R files)
- Reverted outlier removal to old method, but made kept it iterative and in the scree
- Changed binning to the method implemented in the scagnsotics package
- Wrote tests for binning, outlier removal, alpha, and scree functions

# cassowaryr 2.0.10
- Update draw functions so they accept a scree and also draw_alphahull allows you to set the alpha value
- Made the default alpha value clearer (alpha = "rahman")

# cassowaryr 2.0.9
- Updated `scree()` so `outlier_rm` returns revised `del`, `weights`, and `alpha`.
- Added `outlier_rm` and `binner` options to all index functions and `alpha` option to alpha-hull–based indices.
- Added Tina’s name and ORCID to the package description.
- Modified `sc_skinny` and `sc_stringy` to return 1 for perfectly straight-line conditions.


# cassowaryr 2.0.8
- Removed euclid parameter from calculation functions


# cassowaryr 2.0.7
- Rename striated2 to grid and included epsilon parameter to control grid noise tolerance

# cassowaryr 2.0.6
- Change alpha_omega to use mst_weights instead of weights (in the previous version I incorrectly used weights instead of mst_weights)
- Add iterative outlier removal in scree()


# cassowaryr 2.0.5
- Changed `Depends` back to `R (>= 4.0.0)`.
- Removed the Ubuntu oldrel R-CMD-check job from GitHub Actions, since it could not satisfy the previous R (>= 4.0.0) requirement or the new R (>= 4.5.0) dependency and was failing due to dependency errors.


# cassowaryr 2.0.4
- Changed `Depends` to `R (>= 4.5.0)`.
- Implemented hexagonal binning.
- Added hexbin package to Suggests (Description).

# cassowaryr 2.0.3
- Added configurable alpha selection in `scree()`, with support for `"rahman"` (default), `"q90"`, `"omega"`, user-specified numeric values, and user-defined alpha functions. 
- Temporarily skipped `test_alphahull.R` while alpha selection behaviour is being updated (issue).



# cassowaryr 2.0.2 (Current)
- Updated the required version of interp to 1.1-6
- Changed the duplicate check inside the scree calculation to match the check done by interp

# cassowaryr 2.0.1
- Added hex sticker
- Updated the required version of interp to 1.1-4

# cassowaryr 2.0.0

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
