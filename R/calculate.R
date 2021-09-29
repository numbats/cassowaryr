#' Compute scagnostics on all possible scatter plots for the given data
#'
#' @seealso calc_scags
#' @examples
#' # Calculate selected scagnostics
#' calc_scags_wide(pk[,2:5], scags=c("outlying","clumpy","monotonic"))
#'
#' @importFrom magrittr %>%
#' @importFrom progress progress_bar
#' @export
calc_scags_wide <- function(all_data, scags=c("outlying","stringy", "striated", "striated_adjusted", "clumpy", "clumpy_adjusted", "sparse", "skewed", "convex","skinny","monotonic", "splines","dcor"), euclid = TRUE){

  # Check if variables are non-constant
  std_dev <- all_data %>% summarise_all(sd, na.rm=TRUE)
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
    dplyr::group_by(Var1, Var2)%>%
    dplyr::summarise(intermediate_scags(vars=c(Var1, Var2),
                                        data=all_data,
                                        scags=scags, pb))

}

intermediate_scags <- function(vars, data, scags, pb){
  pb$tick()
  x <- dplyr::pull(data, var=vars[[1]])
  y <- dplyr::pull(data, var=vars[[2]])
  return(calc_scags(x,y,scags))
}

#' Compute selected scagnostics on subsets
#'
#' @seealso sc_pairwise
#' @examples
#' # Calculate selected scagnostics on a single pair
#' calc_scags(anscombe$x1, anscombe$y1,
#'   scags=c("monotonic", "outlying","convex"))
#'
#' # Calculate all large number of scagnostics together
#' calc_scags(anscombe$x1, anscombe$y1)
#'
#' # Compute on long form data, or subsets
#' # defined by a categorical variable
#' require(dplyr)
#' datasaurus_dozen %>%
#'   group_by(dataset) %>%
#'   summarise(calc_scags(x,y,
#'     scags=c("monotonic", "outlying", "convex")))
#'
#' @export
calc_scags <- function(x, y, scags=c("outlying", "stringy", "striated", "striated_adjusted", "clumpy", "clumpy_adjusted", "sparse", "skewed", "convex", "skinny", "monotonic", "splines", "dcor")){
  #set all scagnostics to null
  outlying = NULL
  stringy = NULL
  striated = NULL
  striated_adjusted = NULL
  clumpy = NULL
  clumpy_adjusted = NULL
  sparse = NULL
  skewed = NULL
  convex = NULL
  skinny = NULL
  monotonic = NULL
  splines = NULL
  dcor = NULL

  #make original and outlying adjusted scree+mst
  sm_list <- original_and_robust(x,y)

  #scree and mst without outlier removal
  sc_orig <- sm_list$scree_ori
  mst_orig <- sm_list$mst_ori


  #scree and mst with outlier removal
  sc <- sm_list$scree_rob
  mst <- sm_list$mst_rob

  #CALCULATE MST MEASURES
  if("stringy" %in% scags){
    stringy <- sc_stringy.igraph(mst)
  }
  if("striated" %in% scags){
    striated <- sc_striated.igraph(mst, sc)
  }
  if("striated_adjusted" %in% scags){
    striated_adjusted <- sc_striated_adjusted.igraph(mst, sc)
  }
  if("clumpy" %in% scags){
    clumpy <- sc_clumpy.igraph(mst, sc)
  }
  if("clumpy_adjusted" %in% scags){
    clumpy_adjusted <- sc_clumpy_adjusted.igraph(mst, sc)
  }
  if("sparse" %in% scags){
    sparse <- sc_sparse.igraph(mst, sc)
  }
  if("skewed" %in% scags){
    skewed <- sc_skewed.igraph(mst, sc)

  }
  if("outlying" %in% scags){
    outlying <- sc_outlying.igraph(mst_orig, sc_orig)
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
  scagnostic_calcs <- dplyr::tibble("outlying"=outlying,
                             "stringy"=stringy,
                             "striated"=striated,
                             "striated_adjusted" = striated_adjusted,
                             "clumpy"=clumpy,
                             "clumpy_adjusted"=clumpy_adjusted,
                             "sparse"=sparse,
                             "skewed"=skewed,
                             "convex"=convex,
                             "skinny"=skinny,
                             "monotonic"=monotonic,
                             "splines"=splines,
                             "dcor"=dcor)
  return(scagnostic_calcs)
}

