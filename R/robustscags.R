
#function
clumpy_robust <- function(){

}

#function used inside robust clumpy
inner_clumpy <- function(mymst, x){
  #pretty similar to original clumpy with with enough changes that I'm just making a new function
  #input: mst and scree
  #output: Clumpy value and the respective EK and EJ edges

  #lower triangular matrix
  mstmat <- twomstmat(mymst,x)$lowertri

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
  c(clumpy[ind], ej[ind]) #return c(clumpy val, ej weight)

}

#check is need to go deeper on hierarchy bit
split_check <- function(mst, scr){

  #get w value (basically straight from outlying_identify, make own function maybe)
  matlist <- twomstmat(mst,scr)
  mstmat_lt <- matlist$lowertri
  edges <- mstmat_lt[which(mstmat_lt>0)]
  w <- psi(edges)

  #get clumpy value
  clump_ej <- inner_clumpy(mst, scr)

  #no split
  if(clump_ej[2]<w){
    val <- clump_ej[1]
  }
  #split
  if(clump_ej[2]>=w){
    print("big")
    print(clump_ej[2])
    print(w)
  }

}

split_clusters <- function(mst, scr){

}

