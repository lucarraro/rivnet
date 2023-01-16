# extract river with given thrA value
extract_river <- function(river, thrA=NULL, maxReachLength=Inf){

  if (is.null(thrA)){
    thrA <- 0.002*max(river$FD$A,na.rm=T)
  }

  cat("Connectivity at RN level... \n")
  indexRNNodes <- which(river$FD$A >=thrA)

  X_RN <- river$FD$X[indexRNNodes]
  Y_RN <- river$FD$Y[indexRNNodes]
  A_RN <- river$FD$A[indexRNNodes]
  Z_RN <- river$FD$Z[indexRNNodes]
  Length_RN <- river$FD$leng[indexRNNodes]

  ## W, downNode at RN level

  nNodes_RN <- length(indexRNNodes)
  W_RN <- spam(0,nNodes_RN,nNodes_RN)
  ind <- matrix(0,nNodes_RN,2)
  downNode_RN <- numeric(nNodes_RN)
  Slope_RN <-  numeric(nNodes_RN)

  k <- 1
  for (i in 1:nNodes_RN){
    mov <- neigh(river$FD$flowDir[indexRNNodes[i]])
    d <- which(indexRNNodes==(indexRNNodes[i]+mov[1]+mov[2]*river$dimX)) # indices from top-left corner to the right, then next row...
    if (length(d)!=0){
      ind[k, ] <- c(i,d)
      k <- k + 1
      Slope_RN[i] <- (Z_RN[i]-Z_RN[d])/Length_RN[i]
    }
  }
  ind <- ind[-k, ]
  downNode_RN[ind[,1]] <- ind[,2]
  W_RN[ind] <- 1
  rm(ind)
  Outlet_RN <- which(downNode_RN==0)

  # find AG nodes ####
  cat("Connectivity at AG level... \n")
  DegreeIn <- colSums(W_RN)
  DegreeOut <- rowSums(W_RN)
  Confluence <- DegreeIn>1
  Source <- DegreeIn==0
  SourceOrConfluence <- Source|Confluence
  ConfluenceNotOutlet <- Confluence&(downNode_RN!=0)
  ChannelHeads <- SourceOrConfluence  #Source|ConfluenceNotOutlet

  OutletNotChannelHead <- (downNode_RN==0)&(!ChannelHeads)
  IsNodeAG <- SourceOrConfluence|OutletNotChannelHead
  whichNodeAG <- which(IsNodeAG)

  nNodes_AG <- sum(IsNodeAG)
  Length_AG <- numeric(nNodes_AG)
  RN_to_AG <- numeric(nNodes_RN)
  reachID <- 1
  X_AG <- NaN*numeric(nNodes_AG)
  Y_AG <- NaN*numeric(nNodes_AG)
  Z_AG <- NaN*numeric(nNodes_AG)
  A_AG <- NaN*numeric(nNodes_AG)
  while (length(whichNodeAG) != 0){ # explore all AG Nodes
    i <- whichNodeAG[1] # select the first
    RN_to_AG[i] <- reachID
    j <- downNode_RN[i]
    X_AG[reachID] <- X_RN[i]
    Y_AG[reachID] <- Y_RN[i]
    Z_AG[reachID] <- Z_RN[i]
    A_AG[reachID] <- A_RN[i]
    Length_AG[reachID] <- Length_RN[i]
    tmp_length <- Length_RN[i]
    tmp <- NULL
    j0 <- j
    while (!IsNodeAG[j] && j!=0) {
      tmp <- c(tmp, j)
      tmp_length <-  tmp_length + Length_RN[j]
      j_old <- j
      j <- downNode_RN[j]}

    if (tmp_length > maxReachLength){
      n_splits <- ceiling(tmp_length/maxReachLength)
      new_maxLength <- tmp_length/n_splits
      j <- j0
      while (!IsNodeAG[j] && j!=0 && Length_AG[reachID] <= new_maxLength) {
        RN_to_AG[j] <- reachID
        Length_AG[reachID] <-  Length_AG[reachID] + Length_RN[j]
        j_old <- j
        j <- downNode_RN[j]}
      if (Length_AG[reachID] > new_maxLength){
        j <- j_old
        Length_AG[reachID] <-  Length_AG[reachID] - Length_RN[j]
        ChannelHeads[j] <- 1
        whichNodeAG <- c(whichNodeAG,j)}

    } else {
      RN_to_AG[tmp] <- reachID
      Length_AG[reachID] <- tmp_length
    }

    reachID <- reachID + 1
    whichNodeAG <- whichNodeAG[-1]
  }
  nNodes_AG <- length(X_AG)

  # W, downNode at AG level ####

  downNode_AG <- numeric(nNodes_AG)
  W_AG <- spam(0,nNodes_AG,nNodes_AG)
  ind <- matrix(0,nNodes_AG,2)
  reachID <- sum(ChannelHeads) + 1
  for (i in 1:nNodes_RN){
    if (downNode_RN[i] != 0 && RN_to_AG[downNode_RN[i]] != RN_to_AG[i]) {
      downNode_AG[RN_to_AG[i]] <- RN_to_AG[downNode_RN[i]]
      ind[RN_to_AG[i],] <- c(RN_to_AG[i],downNode_AG[RN_to_AG[i]])
    }
  }
  ind <- ind[-which(ind[,1]==0),]
  W_AG[ind] <- 1
  Outlet_AG <- RN_to_AG[Outlet_RN]

  AG_to_RN <- vector("list", nNodes_AG)
  for(i in 1:nNodes_AG) { # attribute river network pixels to fields of the AG_to_FD list
    AG_to_RN[[i]] <- which(RN_to_AG==i)
  }

  FD_to_SC <- NA*numeric(length(flowDir))
  SC_to_FD <- vector("list",nNodes_AG)
  FD_to_SC[indexRNNodes] <- RN_to_AG
  for (i in 1:nNodes_AG){
    SC_to_FD[[i]] <- indexRNNodes[which(RN_to_AG==i)]
  }


  # find FD_to_SC
  drainageArea <- A
  indexFDNodes <- which(drainageArea>0)
  ind_head <- indexFDNodes[which(drainageArea[indexFDNodes]==cellsize^2)]

  for (i in 1:length(ind_head)){
    d <- ind_head[i]
    k <- NA; d_new <- d; sub_d <- numeric(0)
    while (is.na(k)){
      k <- FD_to_SC[d_new]
      if (is.na(k)){
        sub_d <- c(sub_d, d_new)
        mov <- neigh(flowDir[d_new])
        d_new <- d_new + mov[1] + mov[2]*river$dimX # indices from top-left corner to the right, then next row...
      }
    }
    FD_to_SC[sub_d] <- k
    SC_to_FD[[k]] <- c(SC_to_FD[[k]], sub_d)
  }


  # Upstream_RN : list containing IDs of all reaches upstream of each reach (plus reach itself)
  Upstream_RN <- vector("list",nNodes_RN)
  nUpstream_RN <- numeric(nNodes_RN)
  for (i in 1:nNodes_RN){
    UpOneLevel <- which(downNode_RN==i) # find reaches at one level upstream
    Upstream_RN[[i]] <- UpOneLevel      # add them to the list
    while (length(UpOneLevel)!=0) { # continue until there are no more reaches upstream
      ContinuePath <- UpOneLevel # jump 1 level above
      UpOneLevel <- which(downNode_RN %in% ContinuePath) # find reaches at one level upstream
      Upstream_RN[[i]] <- c(Upstream_RN[[i]],UpOneLevel) # add them to the list
    }
    Upstream_RN[[i]] <- c(Upstream_RN[[i]],i)
    nUpstream_RN[i] <- length(Upstream_RN[[i]])
    if ((i %% 100)==0){
      message(sprintf("%.2f%% done\r",100*i/nNodes_RN),appendLF = F)
    }
  }
  message('100.00% done \n',appendLF = F)

  # Upstream_AG : list containing IDs of all reaches upstream of each reach (plus reach itself)
  Upstream_AG <- vector("list",nNodes_AG)
  nUpstream_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    UpOneLevel <- which(downNode_AG==i) # find reaches at one level upstream
    Upstream_AG[[i]] <- UpOneLevel      # add them to the list
    while (length(UpOneLevel)!=0) { # continue until there are no more reaches upstream
      ContinuePath <- UpOneLevel # jump 1 level above
      UpOneLevel <- which(downNode_AG %in% ContinuePath) # find reaches at one level upstream
      Upstream_AG[[i]] <- c(Upstream_AG[[i]],UpOneLevel) # add them to the list
    }
    Upstream_AG[[i]] <- c(Upstream_AG[[i]],i)
    nUpstream_AG[i] <- length(Upstream_AG[[i]])
  }


  # calculate Strahler stream order
  StreamOrder_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    j <- order(nUpstream_AG)[i] # index that explores reaches in a downstream direction
    tmp <- which(downNode_AG==j) # set of reaches draining into j
    if (length(tmp)>0){
      IncreaseOrder <- sum(StreamOrder_AG[tmp]==max(StreamOrder_AG[tmp])) # check whether tmp reaches have the same stream order
      if (IncreaseOrder > 1) {
        StreamOrder_AG[j] <- 1 + max(StreamOrder_AG[tmp]) # if so, increase stream order
      } else {StreamOrder_AG[j] <- max(StreamOrder_AG[tmp])} # otherwise, keep previous stream order
    } else {StreamOrder_AG[j] <- 1} # if j is an headwater, impose StreamOrder = 1
  }


  #print('Length and slope at AG level...',quote=FALSE)
  # Calculate length and slopes of reaches
  #Length_AG <- rep(0,Nnodes_AG)
  Slope_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    #Length_AG[i] <- sum(OCN$FD$leng[AG_to_FD[[i]]])
    Slope_AG[i] <- (Slope_RN[RN_to_AG==i] %*% Length_RN[RN_to_AG==i])/Length_AG[i] # scalar product between vector of slopes and lengths of nodes at RN level belonging to reach i
  }


  nNodes_SC <- nNodes_AG
  Z_SC <- numeric(nNodes_SC)
  Alocal_SC <- numeric(nNodes_SC)
  for (i in 1:nNodes_SC) {
    Z_SC[i] <- mean(river$FD$Z[which(FD_to_SC==i)])
    Alocal_SC[i] <- sum(FD_to_SC==i,na.rm=T)*cellsize^2
  }

  Areach_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG) {
    Areach_AG[i] <- sum(Alocal_SC[Upstream_AG[[i]]])
  }

  # coordinates of AG nodes considered at the downstream end of the respective edge
  XReach <- numeric(nNodes_AG)
  YReach <- numeric(nNodes_AG)
  ZReach <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    tmp <- AG_to_RN[[i]]
    ind <- which(A_RN[tmp]==max(A_RN[tmp]))
    node <- tmp[ind]
    XReach[i] <- X_RN[node]
    YReach[i] <- Y_RN[node]
    ZReach[i] <- Z_RN[node]
  }
  XReach[Outlet_AG] <- NaN
  YReach[Outlet_AG] <- NaN
  ZReach[Outlet_AG] <- NaN

  river$FD[["nNodes"]] <- length(river$FD$X)
  river$FD[["toSC"]] <- FD_to_SC
  river$FD[["toCM"]] <- 1+numeric(length(FD_to_SC))

  river$RN[["A"]] <- A_RN
  river$RN[["downNode"]] <- downNode_RN
  river$RN[["outlet"]] <- which(downNode_RN==0)
  river$RN[["leng"]] <- Length_RN
  river$RN[["upstream"]] <- Upstream_RN
  river$RN[["nUpstream"]] <- nUpstream_RN
  river$RN[["X"]] <- X_RN
  river$RN[["Y"]] <- Y_RN
  river$RN[["toAGReach"]] <- RN_to_AG
  river$RN[["nNodes"]] <- nNodes_RN

  river$AG[["A"]] <- A_AG
  river$AG[["AReach"]] <- Areach_AG
  river$AG[["downNode"]] <- downNode_AG
  river$AG[["upstream"]] <- Upstream_AG
  river$AG[["nUpstream"]] <- nUpstream_AG
  river$AG[["leng"]] <- Length_AG
  river$AG[["outlet"]] <- Outlet_AG
  river$AG[["slope"]] <- Slope_AG
  river$AG[["streamOrder"]] <- StreamOrder_AG
  river$AG[["Upstream"]] <- Upstream_AG
  river$AG[["X"]] <- X_AG
  river$AG[["Y"]] <- Y_AG
  river$AG[["Z"]] <- Z_AG
  river$AG[["W"]] <- W_AG
  river$AG[["nNodes"]] <- nNodes_AG
  pl <- initial_permutation(downNode_AG)
  river$AG[["perm"]] <- pl$perm

  river$SC[["toFD"]] <- SC_to_FD
  river$SC[["A"]] <- Alocal_SC

  river[["thrA"]] <- thrA
  river[["nOutlet"]] <- 1

  cat(sprintf("Total catchment area: %d cells  -   %.2f km2   -  nNodes: %d   -  nLinks: %d \n",
              max(A_RN)/cellsize^2,river$CM$A/1e6,river$RN$nNodes,river$AG$nNodes))

  invisible(river)
}

neigh <- function(dir) {
  mov <- c(0,0)
  switch(dir,
         {mov[1] <- 1; mov[2] <- 0},   # case 1 (E)
         {mov[1] <- 1; mov[2] <- -1},   # case 2 (NE)
         {mov[1] <- 0; mov[2] <- -1},   # case 3 (N)
         {mov[1] <- -1; mov[2] <- -1},  # case 4 (NW)
         {mov[1] <- -1; mov[2] <- 0},  # case 5 (W)
         {mov[1] <- -1; mov[2] <- 1}, # case 6 (SW)
         {mov[1] <- 0; mov[2] <- 1},  # case 7 (S)
         {mov[1] <- 1; mov[2] <- 1})  # case 8 (SE)
  return(mov)
}

# inherited from OCNet
initial_permutation <- function(DownNode){

  Outlet <- which(DownNode==0)
  NodesToExplore <- Outlet # start from outlets
  reverse_perm <- numeric(length(DownNode)) # build permutation vector from outlets to headwaters, then flip it

  k <- 0
  while (length(NodesToExplore)>0){ # continue until all the network has been explored
    k <- k + 1
    node <- NodesToExplore[1] # explore a node
    reverse_perm[k] <- node # assign position in the permutation vector
    NodesToExplore <- NodesToExplore[-1] # remove explored node
    UpNodes <- which(DownNode==node) # find nodes upstream of node
    while (length(UpNodes)>0){ # continue upstream until a headwater is found
      k <- k + 1
      node <- UpNodes[1] # explore first upstream node
      reverse_perm[k] <- node
      if (length(UpNodes)>1){ # if there is a bifurcation upstream, add the other upstream connections at the top of NodesToExplore
        NodesToExplore <- c(UpNodes[2:length(UpNodes)],NodesToExplore)
      }
      UpNodes <- which(DownNode==node)
    }
  }

  perm <- reverse_perm[length(DownNode):1] # flip permutation

  OutList = list(perm=perm,noDAG=0)

  invisible(OutList)
}
