#' Compute robust clumpy scagnostic measure using MST
#'
#' @param x numeric vector of x values
#' @param y numeric vector of y values
#' @return A "numeric" object that gives the plot's robust clumpy score.
#'
#' @examples
#'   require(ggplot2)
#'   require(dplyr)
#'   ggplot(features, aes(x=x, y=y)) +
#'      geom_point() +
#'      facet_wrap(~feature, ncol = 5, scales = "free")
#'   features %>% group_by(feature) %>% summarise(clumpy = sc_clumpy_r(x,y))
#'   sc_clumpy_r(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'
#' @export
sc_clumpy_r <- function(x, y) UseMethod("sc_clumpy_r")

#' @rdname sc_clumpy_r
#' @export
sc_clumpy_r.default <- function(x, y){
  sc <- scree(x, y)
  sc_clumpy_r.scree(sc)
}

#' @rdname sc_clumpy_r
#' @export
sc_clumpy_r.scree <- function(x, y=NULL) {
  #generate vector of MST edges
  mst <- gen_mst(x$del, x$weights)
  sc_clumpy_r.igraph(mst, x)
}

#' @rdname sc_clumpy_r
#' @export
sc_clumpy_r.igraph <- function(x, y){
  mst_lt <- twomstmat(x,y)$lowertri
  vals <- outside_cluster(mst_lt)
  n <- length(which(mst_lt>0))
  sum(vals)/n
}

# function used inside robust clumpy_r
inner_clumpy <- function(mstmat){
  #pretty similar to original clumpy with with enough changes that I'm just making a new function
  #input: mst (lower triangular matrix) and scree
  #output: Clumpy value and the respective EK and EJ edges

  #make index variables to iterate through
  matind <- which(mstmat>0)
  rows <- matind %% length(mstmat[,1])
  cols <- (matind-1) %/% length(mstmat[,1]) +1
  clumpy <- rep(0,length(matind))
  ej <- rep(0,length(matind))

  for(j in seq(length(rows))){
    #set mst we are going to remove all the edges from
    mst_ej <- mstmat

    #have two clusters sprawling out from the deleted edge (inex i)
    c1rowcol <- rows[j]
    c2rowcol <- cols[j]

    #get weight of ej
    ej_weight <- mst_ej[matind[j]]

    #remove ej and all values in mst that are greater than ej
    mst_ej[which(mst_ej>=ej_weight)] = 0

    #remake index objects to edit within iteration
    matind_ej <- which(mst_ej>0)
    rows_ej <- matind_ej %% length(mst_ej[,1])
    cols_ej <- (matind_ej-1) %/% length(mst_ej[,1]) +1

    #initialise variable that checks if clusters have changed
    whilecheck <- c(c1rowcol,c2rowcol)
    quitloop = 0

    #add in new indices until both clusters are full
    while(quitloop==0){

      #find matches in rows/columns to join connected vertices
      c1rowcol <- unique(c(c1rowcol,
                           rows_ej[which(cols_ej%in%c1rowcol)],
                           cols_ej[which(rows_ej%in%c1rowcol)]))
      c2rowcol <- unique(c(c2rowcol,
                           rows_ej[which(cols_ej%in%c2rowcol)],
                           cols_ej[which(rows_ej%in%c2rowcol)]))

      #check if indices are done updating
      if(setequal(whilecheck, c(c1rowcol,c2rowcol))){
        quitloop=1
      }

      #update while loop check
      whilecheck = c(c1rowcol,c2rowcol)

    }

    #get matrix indices of each of the clusters
    c1ind <- unique(c(matind_ej[which(cols_ej%in%c1rowcol)],
                      matind_ej[which(rows_ej%in%c1rowcol)]))
    c2ind <- unique(c(matind_ej[which(cols_ej%in%c2rowcol)],
                      matind_ej[which(rows_ej%in%c2rowcol)]))

    #get edge weights for each cluster
    c1weights <- mst_ej[c1ind]
    c2weights <- mst_ej[c2ind]

    #set K and M weights
    ek_weight <- max(c1weights, na.rm=TRUE) #max(c(0,c1weights))
    em_weight <- max(c2weights, na.rm=TRUE) #max(c(0,c2weights))

    #calculate this clumpy value
    numer1 <- length(c1ind)*ek_weight
    numer2 <- length(c2ind)*em_weight
    denom <- (length(c1ind)+length(c2ind))*ej_weight

    #maintain vectors
    clumpy[j] <- 1 - (numer1+numer2)/denom #ifelse(denom==0, 0, 1 - (numer1+numer2)/denom) #set =0 if denom=0
    ej[j] <- ej_weight

  }
  # return final clumpy measure
  ind <- which(clumpy==max(clumpy, na.rm=TRUE)) #which.max(clumpy)
  c(clumpy[ind], ej[ind], matind[ind]) #return c(clumpy val, ej weight, ind in matrix)
}

#check is need to go deeper on hierarchy bit
outside_cluster <- function(mstmat_lt){

  # find outlier threshold
  edges <- mstmat_lt[which(mstmat_lt>0)] #get edges
  w <- psi(edges) #set threshold for outliers

  # TODO: REMOVE OUTLIERS


  #get clumpy value
  clump_ej <- inner_clumpy(mstmat_lt)

  #check if leaf actualll has clumpy value
  val=0

  if(length(clump_ej)>0){
    #no split
    if(clump_ej[2]<w){
      val <- length(edges)*clump_ej[1]
    }
    #split
    if(clump_ej[2]>=w){
      #split into two MST
      newclusters <- split_clusters(mstmat_lt, clump_ej[3])
      val <- c(outside_cluster(newclusters$mat1), outside_cluster(newclusters$mat2))
    }
  }

  val #retun value
}

split_clusters <- function(mst_lt, ind, splitforclumpy=FALSE){
  #make index variables to iterate through
  matind <- which(mst_lt>0)
  j <- which(matind == ind)

  #have two clusters sprawling out from the deleted edge (inex i)
  c1rowcol <- (matind %% length(mst_lt[,1]))[j]
  c2rowcol <- ((matind-1) %/% length(mst_lt[,1]) +1)[j]

  #remove ej
  if (splitforclumpy==FALSE){
    mst_lt[ind] = 0
  }
    #and all values in mst that are greater than ej if inside clumpy calc
  if (splitforclumpy==TRUE){
    ej_weight <- mst_lt[ind]
    mst_lt[which(mst_lt >= ej_weight)] = 0
  }

  #remake index objects to edit within iteration
  matind_ej <- which(mst_lt>0)
  rows_ej <- matind_ej %% length(mst_lt[,1])
  cols_ej <- (matind_ej-1) %/% length(mst_lt[,1]) +1

  #initialise variable that checks if clusters have changed
  whilecheck <- c(c1rowcol,c2rowcol)
  quitloop = 0

  while(quitloop==0){

    # find matches in rows/columns to join connected vertices
    c1rowcol <- unique(c(c1rowcol,
                         rows_ej[which(cols_ej%in%c1rowcol)],
                         cols_ej[which(rows_ej%in%c1rowcol)]))
    c2rowcol <- unique(c(c2rowcol,
                         rows_ej[which(cols_ej%in%c2rowcol)],
                         cols_ej[which(rows_ej%in%c2rowcol)]))

    # check if indices are done updating
    if(setequal(whilecheck, c(c1rowcol,c2rowcol))){
      quitloop=1
    }

    #update while loop check
    whilecheck = c(c1rowcol,c2rowcol)

  }

  #get matrix indices of each of the clusters
  c1ind <- unique(c(matind_ej[which(cols_ej%in%c1rowcol)],
                    matind_ej[which(rows_ej%in%c1rowcol)]))
  c2ind <- unique(c(matind_ej[which(cols_ej%in%c2rowcol)],
                    matind_ej[which(rows_ej%in%c2rowcol)]))
  #make two different adjacency matrices based on two clusters
  mst1 <- mst_lt
  mst1[c2ind] = 0
  mst2 <- mst_lt
  mst2[c1ind] = 0
  list(mat1 = mst1,
       mat2 = mst2)
}

