river_to_shapefile <- function(river, filename, atts=NULL,
                               EPSG=NULL, ...){

  dots <- list(...)
  Id <- X <- Y  <- NULL

  for (i in 1:river$AG$nNodes){
    j <- which(river$RN$toAG==i)
    path <- j
    while (river$RN$toAG[j]==i | river$RN$toAG[j]==0){
      j <- river$RN$downNode[j]
      if (j==0) {path <- c(path, path); break} # if we are at the outlet
      path <- c(path, j)
    }
    X <- c(X, river$RN$X[path])
    Y <- c(Y, river$RN$Y[path])
    Id <- c(Id, rep(i,length(path)))

  }
  if (!is.null(atts)){
    atts_df <- data.frame(matrix(0, river$AG$nNodes, length(atts)))
    colnames(atts_df) <- atts
    for (ff in atts){
      atts_df[[ff]] <- river$AG[[ff]]
    }
  } else {atts_df <- NULL}

  mm <- matrix(0, length(X), 5)
  mm[,1] <- Id; mm[,2] <- 1
  mm[,3] <- X;  mm[,4] <- Y
  colnames(mm) <- c("object","part","X","Y","hole")

  v <- vect(mm, "lines", atts=atts_df)
  if (!is.null(EPSG)){
    str <- paste0("epsg:",EPSG)
    crs(v) <- str
  }
  dots$x <- v
  dots$filename <- filename

  do.call(writeVector, dots)

}

