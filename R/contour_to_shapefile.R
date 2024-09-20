contour_to_shapefile <- function(river, filename,
                                 EPSG=NULL, ...){

  dots <- list(...)
  k <- 0
  Id <- X <- Y <- Id2 <- NULL
  n1 <- length(river$CM$XContour)
  for (i1 in 1:n1){
    n2 <- length(river$CM$XContour[[i1]])
    for (i2 in 1:n2){
      ll <- length(river$CM$XContour[[i1]][[i2]])

      X <- c(X, river$CM$XContour[[i1]][[i2]])
      Y <- c(Y, river$CM$YContour[[i1]][[i2]])
      Id <- c(Id, rep(i1, ll))
      Id2 <- c(Id2, rep(i2, ll))
    }
  }

  mm <- matrix(0, length(X), 5)
  mm[,1] <- Id; mm[,2] <- Id2
  mm[,3] <- X;  mm[,4] <- Y
  colnames(mm) <- c("object","part","x","y","hole")

  v <- vect(mm, "polygons")
  if (!is.null(EPSG)){
    str <- paste0("epsg:",EPSG)
    crs(v) <- str
  }
  dots$x <- v
  dots$filename <- filename

  do.call(writeVector, dots)

}
