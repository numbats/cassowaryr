

inside_cluster <- function(mst, sc){
  #mst matrix
  matlist <- twomstmat(mst,scr)#function for matrix version of mst
  mstmat <- matlist$mat #mst matrix
  mstmat_lt <- matlist$lowertri #mst lower triangular matrix

  #calculate w
  edges <- mstmat_lt[which(mstmat>0)]
  w <- psi(edges)

  #remove outliers
  outliers <- outlying_identify(mst, sc)

  #calculate clumpy values
  value <- clumpy_value()
  #return ej, clumpy val, w
}

clumpy_value <- function(mymst, x){
  #input: mst and scree
  #output: clumpy value
  #lower triangular matrix
  mstmat <- twomstmat(mymst,x)$lowertri

  #get cols and rows of each value
  edges <- which(mstmat>0)
  rows <- edges %% length(mstmat[,1])
  cols <- (edges-1) %/% length(mstmat[,1]) +1
  print(edges)
  print(rows)
  print(cols)
  #clumpy <- rep(0,length(edges))

  #for(i in seq(length(edges))){
  #  #set value for edge in consideration
  #  jval <- mstmat[edges[i]]
    #reset inner loop values and remove target edge
  #  inedges <- edges[-i]
  #  inrows <- rows[-i]
  #  incols <- cols[-i]
    #make cluster with remaining edges

  #  kval <- ifelse(sum(cluster1)<sum(cluster2), pmax(cluster1), pmax(cluster2))
  #  clumpy[i] <- 1- (kval/jval)
  #}
  #clumpy <- clumpy[which(!is.na(clumpy))]
  #max(clumpy)
}

two_clusters <- function(inedges,inrows,incols){
  #set value for edge in consideration
  jval <- mstmat[edges[i]]
  #reset inner loop values and remove target edge
  inedges <- edges[-i]
  inrows <- rows[-i]
  incols <- cols[-i]

  #input: the weights, row index and column index of edges
  #output: the index of cluster 1 in the MST matrix
  for(j in seq(length(inedges))){
    #start cluster 1 with
    if(j==1){
      ind=1
    }
    checkvec <- c(inrows[ind],incols[ind]) #lower triang matrix so need to check rows and cols
    if((inrows[j]%in%checkvec | incols[j]%in%checkvec)& (j!=1)){
      ind <- c(ind,j)
    }
  }
  ind
}

