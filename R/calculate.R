#' Compute scagnostics on all possible scatter plots for the given data
#'
#' It is quite common to have data in a wide format that is not suitable to
#' feed into the calc_scags function that would need a long format. To save
#' users time and energy we also provide a wide version of the calc_scags
#' function. This function will compute all selected scagnostics for every
#' pair wise set of variables in the data frame.
#'
#' @inheritParams scree
#' @param all_data tibble of wide multivariate data
#' @param scags collection of strings matching names of
#' scagnostics to calculate: outlying, stringy, striated,
#' grid, striped, clumpy, clumpy2, sparse, skewed, convex,
#' skinny, monotonic, splines, dcor
#'
#' @return A data frame that gives the scagnostic scores for every possible
#' pair of variables.
#'
#' @seealso calc_scags
#' @examples
#' # Calculate selected scagnostics
#' data(pk)
#' calc_scags_wide(pk[,2:5], scags=c("outlying","monotonic"))
#'
#' @importFrom progress progress_bar
#' @export
calc_scags_wide <- function(all_data, scags=c("outlying", "stringy", "striated",
                                              "clumpy", "sparse", "skewed",
                                              "convex", "skinny", "monotonic"),
                            out.rm= TRUE, binner = "hex", alpha = "rahman"){

  if("striated2" %in% scags){
    warning("Please use grid instead of striated2")
    scags[which(scags=="striated2")] <- "grid"
    }

  # Check for typos/misspellings in scags list
  validscags <- c("outlying", "stringy", "striated", "grid", "clumpy", "clumpy2", "sparse", "skewed", "convex", "skinny", "monotonic", "splines", "dcor")
  validpased <- match.arg(scags, validscags, several.ok=TRUE)
  typo <- scags[!scags %in% validpased]
  if(length(typo)>0){
    warning(paste0("You passed ", typo, " to the scags option. That is not a scagnsotic, Did you make a typo?"))
    scags <- validpased #replace scags list with only valid scags
  }

  # Check if variables are non-constant
  Var1 <- Var2 <- NULL
  std_dev <- all_data |> dplyr::summarise_all(stats::sd, na.rm=TRUE)
  keep <- names(std_dev)[std_dev > 0]
  drop <- names(std_dev)[!(names(std_dev) %in% keep)]
  if (length(drop) > 0) {
    message("These variables have constant values and are removed: \n")
    for (i in 1:length(drop))
      message(drop[i])
    all_data <- all_data[, keep]
  }

  # make a dataset of all pairwise variable combinations
  all_combs <- expand.grid(colnames(all_data), colnames(all_data)) |>
    dplyr::filter(!(Var1==Var2))

  # get rid of reversed duplicates
  all_combs <- all_combs[!duplicated(apply(all_combs, 1, function(x) paste(sort(x), collapse=''))),]

  # set up progress bar (ticks in intermediate scags)
  num_ticks <- length(all_combs$Var1)
  pb <- progress_bar$new(format = "[:bar] (:percent) eta :eta", total = num_ticks)

  # calculate scagnostics
  all_combs |>
    dplyr::group_by(Var1, Var2) |>
    dplyr::summarise(intermediate_scags(vars=c(Var1, Var2),
                                        data=all_data,
                                        scags=scags, out.rm, pb)) |>
    dplyr::ungroup()

}

intermediate_scags <- function(vars, data, scags, out.rm, pb){
  pb$tick()
  x <- dplyr::pull(data, var=vars[[1]])
  y <- dplyr::pull(data, var=vars[[2]])
  return(calc_scags(x, y, scags, out.rm))
}

#' Compute selected scagnostics on subsets
#'
#' This function allows you to calculate a large number of scagnostics quickly
#' and efficiently. While the individual scagnostic calculation functions
#' (sc_*) are good for looking at a single scagnostic, it is inefficient when
#' computing more than one scagnostic. This is becasue the sc_* functions
#' recompute the graph object for each plot and scagnostic pair, even if
#' the graph objects are unchanged. Additionally, typing the functions over
#' and over again quickly becomes tedious. For this reason, we have the
#' calc_scags function that will reuse the same graph
#' object for all scagnostics.
#'
#' @inheritParams scree
#' @param scags collection of strings matching names of
#' scagnostics to calculate: outlying, stringy, striated,
#' grid, striped, clumpy, clumpy2, sparse, skewed, convex,
#' skinny, monotonic, splines, dcor
#'
#' @return A data frame with all selected scagnostic values for a particular x,
#' y pair.
#'
#' @seealso calc_scags_wide
#' @examples
#' # Calculate selected scagnostics on a single pair
#' calc_scags(anscombe$x1, anscombe$y1, scags=c("monotonic", "outlying"))
#'
#' # Compute on long form data, or subsets
#' # defined by a categorical variable
#' require(dplyr)
#' datasaurus_dozen |>
#'   group_by(dataset) |>
#'   summarise(calc_scags(x,y, scags=c("monotonic", "outlying", "convex")))
#'
#' @export
calc_scags <- function(x, y, scags=c("outlying", "stringy", "striated",
                                     "clumpy", "sparse", "skewed",
                                     "convex", "skinny", "monotonic"),
                       out.rm=TRUE, binner = "hex", alpha = "rahman"){
  #set all scagnostics to null
  outlying = NULL
  stringy = NULL
  striated = NULL
  grid = NULL
  clumpy = NULL
  clumpy2 = NULL
  sparse = NULL
  skewed = NULL
  convex = NULL
  skinny = NULL
  monotonic = NULL
  splines = NULL
  dcor = NULL
  striped = NULL
  # replace striated2 with grid
  if("striated2" %in% scags){
    warning("'striated2' is no longer an available scagnostic, please use 'grid' instead")
    scags[which(scags=="striated2")] <- "grid"
    }

  # get list for each object type for easier computation
  xy_scags <- c("monotonic", "splines", "dcor", "striped")
  ahull_scags <- c("convex", "skinny")
  mst_scags <- c("outlying", "stringy", "stringy2", "striated", "grid",
                "clumpy", "clumpy2", "sparse", "skewed")

  # Check for typos/misspellings in scags list
  validscags <- c(mst_scags, ahull_scags, xy_scags)
  validpased <- match.arg(scags, validscags, several.ok=TRUE)
  typo <- scags[!scags %in% validpased]
  if(length(typo)>0){
    warning(paste0("You passed ", typo, " to the scags option. That is not a scagnsotic, Did you make a typo?"))
    scags <- validpased #replace scags list with only valid scags
  }

  # 1) Calculate association measures using x,y values
  if("monotonic" %in% scags){
    monotonic <- sc_monotonic(x,y)
  }
  if("splines" %in% scags){
    splines <- sc_splines(x,y)
  }
  if("dcor" %in% scags){
    dcor <- sc_dcor(x,y)
  }
  # Striped, special index
  if("striped" %in% scags){
    striped <- sc_striped(x,y)
  }

  # 2) Calculate outlying measure by itself
  # (it can be calculated with other MST scags if out.rm = FALSE)
  if(out.rm == TRUE & "outlying" %in% scags){
    outlying <- sc_outlying(x,y, binner = binner)

    # remove from scag list so it is not recomputed
    scags <- setdiff(scags, "outlying")
  }

  # 2) Implement delauney triangulation if we have graph based scags
  if(length(scags) == 0){
    break
  } else{
    scr <- scree(x, y, out.rm = out.rm, binner = binner, alpha = alpha)
  }

  # 3) If we have MST scags, calculate MST
  if(any(scags %in% mst_scags)){
    mst <- gen_mst(scr$del, scr$weights)

    # 4) Then calculate MST scagnostics
    if("outlying" %in% scags){
      stringy <- sc_stringy(mst)
    }
    if("stringy" %in% scags){
      stringy <- sc_stringy(mst)
    }
    if("striated" %in% scags){
      striated <- sc_striated(mst, scr)
    }
    if("grid" %in% scags){
      grid <- sc_grid(mst, scr)
    }
    if("clumpy" %in% scags){
      clumpy <- sc_clumpy(mst)
    }
    if("clumpy2" %in% scags){
      clumpy2 <- sc_clumpy2(mst, scr)
    }
    if("sparse" %in% scags){
      sparse <- sc_sparse(mst)
    }
    if("skewed" %in% scags){
      skewed <- sc_skewed(mst)
    }
    # Remove MST scags from list
    scags <- setdiff(scags, mst_scags)
  }

  # 4) If any hull based scags
  if(any(scags %in% ahull_scags)){

    # Calculate alpha and convex hulls
    ahull <- gen_alpha_hull(scr$del, scr$alpha)
    chull <- gen_conv_hull(scr$del)

    # Then calculate the scagnostics
    if("convex" %in% scags){
      convex <- sc_convex.list(chull,ahull)
    }
    if("skinny" %in% scags){
      skinny <- sc_skinny.list(ahull)
    }
  }

  # return table of scagnostic values
  scagnostic_calcs <- dplyr::tibble("outlying"=outlying,
                             "stringy"=stringy,
                             "striated"=striated,
                             "grid" = grid,
                             "clumpy"=clumpy,
                             "clumpy2"=clumpy2,
                             "sparse"=sparse,
                             "skewed"=skewed,
                             "convex"=convex,
                             "skinny"=skinny,
                             "monotonic"=monotonic,
                             "splines"=splines,
                             "dcor"=dcor,
                             "striped"=striped)
  return(scagnostic_calcs)
}

#' Summary computations for scagnostic data
#'
#' These functions suggests a summary statistic that can be found using
#' the scag calculations provided by calc_scags. The top_pair function
#' finds the top pair of variables for each of the scagnostics, while
#' top_scag finds the highest value scagnostic for each pair of variables.
#' While these computations
#' are relatively straight forward for any R user to compute themselves,
#' including these summary function in the package simultaneously streamlines
#' a common calculation made with the scagnostic data and suggests this
#' summary to new users of the package.
#'
#' @param scags_data A dataset of scagnostic values that was returned by
#' calc_scags or calc_scags_wide
#'
#' @return A data frame. For top_pair, each row will represent a scagnostic
#' with its highest pair. For top_scag, each row will represent a pair of
#' variables with its highest valued scagnostic.

#' @examples
#' require(dplyr)
#' # calculate scag data
#' scag_data <- datasaurus_dozen |>
#'   group_by(dataset) |>
#'   summarise(calc_scags(x,y, scags=c("monotonic", "outlying", "convex")))
#'
#' # Calculate top_pair
#' scag_data |>
#'   top_pair()
#'
#' # Calculate top_scag
#' scag_data |>
#'   top_scag()
#'
#'
#' @seealso calc_scags calc_scags_wide
#' @name top_functions

#' @rdname top_functions
#' @export
top_pair <- function(scags_data){
  value <- scag <- NULL
  validscags <- c("outlying", "stringy", "striated", "grid", "clumpy",
                  "clumpy2", "sparse", "skewed", "convex", "skinny",
                  "monotonic", "splines", "dcor", "striped")
  scags_data |>
    tidyr::pivot_longer(tidyselect::any_of(validscags), names_to = "scag",
                        values_to = "value") |>
    dplyr::arrange(dplyr::desc(value)) |>
    dplyr::group_by(scag) |>
    dplyr::slice_head(n=1) |>
    dplyr::ungroup()
}


#' @rdname top_functions
#' @export
top_scag <- function(scags_data){
  value <- scag <- NULL
  validscags <- c("outlying", "stringy", "striated", "grid", "clumpy", "clumpy2", "sparse", "skewed", "convex", "skinny", "monotonic", "splines", "dcor", "striped")
  scags_data |>
    tidyr::pivot_longer(tidyselect::any_of(validscags), names_to = "scag", values_to = "value") |>
    dplyr::arrange(dplyr::desc(value)) |>
    dplyr::group_by(dplyr::across(-c(scag,value))) |>
    dplyr::slice_head(n=1) |>
    dplyr::ungroup()
}


