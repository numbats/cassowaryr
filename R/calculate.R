#' Compute scagnostics on all possible scatter plots for the given data
#'
#' @param all_data tibble of multivariate data on which to compute scagnostics
#' @param scags collection of strings matching names of
#' scagnostics to calculate: outlying, stringy, striated,
#' striated2, striped, clumpy, clumpy2, sparse, skewed, convex,
#' skinny, monotonic, splines, dcor
#' @param euclid logical indicator to use Euclidean distance
#' @param out.rm logical indicator to indicate if outliers should be removed before calculating non outlying measures
#' @return A data frame that gives the data's scagnostic scores for each possible variable combination.
#' @seealso calc_scags
#' @examples
#' # Calculate selected scagnostics
#' data(pk)
#' calc_scags_wide(pk[,2:5], scags=c("outlying","monotonic"))
#'
#' @importFrom magrittr %>%
#' @importFrom progress progress_bar
#' @export
calc_scags_wide <- function(all_data, scags=c("outlying", "stringy", "striated", "striated2", "clumpy", "clumpy2", "sparse", "skewed", "convex", "skinny", "monotonic", "splines", "dcor"), out.rm= TRUE, euclid = FALSE){

  # Check for typos/misspellings in scags list
  validscags <- c("outlying", "stringy", "striated", "striated2", "clumpy", "clumpy2", "sparse", "skewed", "convex", "skinny", "monotonic", "splines", "dcor")
  validpased <- match.arg(scags, validscags, several.ok=TRUE)
  typo <- scags[!scags %in% validpased]
  if(length(typo)>0){
    warning(paste0("You passed ", typo, " to the scags option. That is not a scagnsotic, Did you make a typo?"))
    scags <- validpased #replace scags list with only valid scags
  }

  # Check if variables are non-constant
  Var1 <- Var2 <- NULL
  std_dev <- all_data %>% dplyr::summarise_all(stats::sd, na.rm=TRUE)
  keep <- names(std_dev)[std_dev > 0]
  drop <- names(std_dev)[!(names(std_dev) %in% keep)]
  if (length(drop) > 0) {
    message("These variables have constant values and are removed: \n")
    for (i in 1:length(drop))
      message(drop[i])
    all_data <- all_data[, keep]
  }

  # make a dataset of all pairwise variable combinations
  all_combs <- expand.grid(colnames(all_data), colnames(all_data)) %>%
    dplyr::filter(!(Var1==Var2))

  # get rid of reversed duplicates
  all_combs <- all_combs[!duplicated(apply(all_combs, 1, function(x) paste(sort(x), collapse=''))),]

  # set up progress bar (ticks in intermediate scags)
  num_ticks <- length(all_combs$Var1)
  pb <- progress_bar$new(format = "[:bar] (:percent) eta :eta", total = num_ticks)

  # calculate scagnostics
  all_combs %>%
    dplyr::group_by(Var1, Var2) %>%
    dplyr::summarise(intermediate_scags(vars=c(Var1, Var2),
                                        data=all_data,
                                        scags=scags, out.rm, euclid, pb)) %>%
    dplyr::ungroup()

}

intermediate_scags <- function(vars, data, scags, out.rm, euclid, pb){
  pb$tick()
  x <- dplyr::pull(data, var=vars[[1]])
  y <- dplyr::pull(data, var=vars[[2]])
  return(calc_scags(x, y, scags, out.rm, euclid))
}

#' Compute selected scagnostics on subsets
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param scags collection of strings matching names of
#' scagnostics to calculate: outlying, stringy, striated,
#' striated2, striped, clumpy, clumpy2, sparse, skewed, convex,
#' skinny, monotonic, splines, dcor
#' @param out.rm logical indicator to indicate if outliers should be removed before calculating non outlying measures
#' @param euclid logical indicator to use Euclidean distance
#' @return A data frame that gives the single plot's scagnostic score.
#' @seealso calc_scags_wide
#' @examples
#' # Calculate selected scagnostics on a single pair
#' calc_scags(anscombe$x1, anscombe$y1, scags=c("monotonic", "outlying"))
#'
#' # Compute on long form data, or subsets
#' # defined by a categorical variable
#' require(dplyr)
#' datasaurus_dozen %>%
#'   group_by(dataset) %>%
#'   summarise(calc_scags(x,y, scags=c("monotonic", "outlying", "convex")))
#'
#' @export
calc_scags <- function(x, y, scags=c("outlying", "stringy", "striated", "striated2", "clumpy", "clumpy2", "sparse", "skewed", "convex", "skinny", "monotonic", "splines", "dcor"), out.rm=TRUE, euclid=FALSE){
  #set all scagnostics to null
  outlying = NULL
  stringy = NULL
  striated = NULL
  striated2 = NULL
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
  euclid_dist= NULL

  # Check for typos/misspellings in scags list
  validscags <- c("outlying", "stringy", "striated", "striated2", "clumpy", "clumpy2", "sparse", "skewed", "convex", "skinny", "monotonic", "splines", "dcor", "striped")
  validpased <- match.arg(scags, validscags, several.ok=TRUE)
  typo <- scags[!scags %in% validpased]
  if(length(typo)>0){
    warning(paste0("You passed ", typo, " to the scags option. That is not a scagnsotic, Did you make a typo?"))
    scags <- validpased #replace scags list with only valid scags
  }

  # Remove missing values
  d <- tibble::tibble(x=x, y=y)
  d <- d[stats::complete.cases(d),]
  if (length(x) > nrow(d))
    message("WARNING: ", length(x)-nrow(d), " observations in have been removed. \n")
  x <- d$x
  y <- d$y

  # Make original and outlying adjusted scree+mst
  sm_list <- original_and_robust(x,y)

  if(is.null(sm_list)){
    # this is null when one of the variables is constant after outlier removal
    return(
      dplyr::tibble("outlying"=NA,
                    "stringy"=NA,
                    "striated"=NA,
                    "striated2" = NA,
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

  #scree and mst without outlier removal
  sc_orig <- sm_list$scree_ori
  mst_orig <- sm_list$mst_ori

  #scree and mst with outlier removal
  sc <- sm_list$scree_rob
  mst <- sm_list$mst_rob

  #calculate outlying meausre
  if("outlying" %in% scags){
    outlying <- sc_outlying(mst_orig, sc_orig)
  }

  #check is outliers should be removed
  if(out.rm == FALSE){
    sc <- sc_orig
    mst <- mst_orig
  }

  #CALCULATE MST MEASURES
  if("stringy" %in% scags){
    stringy <- sc_stringy(mst)
  }
  if("striated" %in% scags){
    striated <- sc_striated(mst, sc)
  }
  if("striated2" %in% scags){
    striated2 <- sc_striated2(mst, sc)
  }
  if("clumpy" %in% scags){
    clumpy <- sc_clumpy(mst, sc)
  }
  if("clumpy2" %in% scags){
    clumpy2 <- sc_clumpy2(mst, sc)
  }
  if("sparse" %in% scags){
    sparse <- sc_sparse(mst, sc)
  }
  if("skewed" %in% scags){
    skewed <- sc_skewed(mst, sc)

  }

  #CALCULATE ALPHA HULL MEASURES
  if(any("convex" %in% scags, "skinny" %in% scags)){
    chull <- gen_conv_hull(sc$del)
    ahull <- gen_alpha_hull(sc$del, sc$alpha)
    if("convex" %in% scags){
      convex <- sc_convex.list(chull,ahull)
    }
    if("skinny" %in% scags){
      skinny <- sc_skinny.list(ahull)
    }}

  #CALCULATE ASSOCIATION MEASURES
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

  #calculate euclidean distance
  if(euclid==TRUE){
    scagvect <- c(outlying, stringy, striated, striated2, clumpy, clumpy2,
                  sparse, skewed, convex, skinny, monotonic, splines, dcor, striped)
    euclid_dist <- sqrt(sum(scagvect^2))
  }

  scagnostic_calcs <- dplyr::tibble("outlying"=outlying,
                             "stringy"=stringy,
                             "striated"=striated,
                             "striated2" = striated2,
                             "clumpy"=clumpy,
                             "clumpy2"=clumpy2,
                             "sparse"=sparse,
                             "skewed"=skewed,
                             "convex"=convex,
                             "skinny"=skinny,
                             "monotonic"=monotonic,
                             "splines"=splines,
                             "dcor"=dcor,
                             "striped"=striped,
                             "euclid_dist" = euclid_dist)
  return(scagnostic_calcs)
}

#' Calculate the top pair of variables or group for each scagnostic
#' @param scags_data A dataset of scagnostic values that was returned by calc_scags or calc_scags_wide
#' @return A data frame where each row is a scagnostic with its highest pair and the associated value
#' @examples
#' #an example using calc_scags
#' require(dplyr)
#' datasaurus_dozen %>%
#'   group_by(dataset) %>%
#'   summarise(calc_scags(x,y, scags=c("monotonic", "outlying", "convex"))) %>%
#'   top_scags()
#'  #an example using calc_scags_wide
#'  data(pk)
#'  scags_data <- calc_scags_wide(pk[,2:5], scags=c("outlying","clumpy","monotonic"))
#'  top_scags(scags_data)
#' @seealso calc_scags calc_scags_wide top_pairs
#' @importFrom magrittr %>%
#' @export
top_scags <- function(scags_data){
  value <- scag <- NULL
  validscags <- c("outlying", "stringy", "striated", "striated2", "clumpy", "clumpy2", "sparse", "skewed", "convex", "skinny", "monotonic", "splines", "dcor", "striped")
  scags_data %>%
    tidyr::pivot_longer(tidyselect::any_of(validscags), names_to = "scag", values_to = "value") %>%
    dplyr::arrange(dplyr::desc(value)) %>%
    dplyr::group_by(scag) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::ungroup()
}

#' Calculate the top scagnostic for each pair of variables
#' @param scags_data A dataset of scagnostic values that was returned by calc_scags or calc_scags_wide
#' @return A data frame where each row is a each scatter plot, its highest valued scagnostic, and its respective value
#' @examples
#' #an example using calc_scags
#' require(dplyr)
#' datasaurus_dozen %>%
#'   group_by(dataset) %>%
#'   summarise(calc_scags(x,y, scags=c("monotonic", "outlying", "convex"))) %>%
#'   top_pairs()
#'  #an example using calc_scags_wide
#'  data(pk)
#'  scags_data <- calc_scags_wide(pk[,2:5], scags=c("outlying","clumpy","monotonic"))
#'  top_pairs(scags_data)
#' @seealso calc_scags calc_scags_wide top_scags
#' @importFrom magrittr %>%
#' @export
top_pairs <- function(scags_data){
  value <- scag <- NULL
  validscags <- c("outlying", "stringy", "striated", "striated2", "clumpy", "clumpy2", "sparse", "skewed", "convex", "skinny", "monotonic", "splines", "dcor", "striped")
  scags_data %>%
    tidyr::pivot_longer(tidyselect::any_of(validscags), names_to = "scag", values_to = "value") %>%
    dplyr::arrange(dplyr::desc(value)) %>%
    dplyr::group_by(dplyr::across(-c(scag,value))) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::ungroup()
}


