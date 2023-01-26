locate_site <- function(X,Y,river,
                        euclidean=TRUE,
                        showPlot=FALSE,
                        xlim=NULL,
                        ylim=NULL){

  if (!("RN" %in% names(river))){
    stop('Missing fields in river. You should run aggregate_river prior to locate_site.')
  }

  distVec_FD <- sqrt((river$FD$X-X)^2 + (river$FD$Y-Y)^2)
  indFD <- which(distVec_FD == min(distVec_FD))[1]

  Xgrid <- river$FD$X[indFD]
  Ygrid <- river$FD$Y[indFD]

  if (euclidean){
    # find closest site as the crow flies
    distVec <- sqrt((river$RN$X-Xgrid)^2 + (river$RN$Y-Ygrid)^2)
    distanz <- min(distVec)
    RNnode <- which(distVec == distanz)[1]
  } else {
    distanz <- 0
    # follow downstream direction
    j <- indFD
    while (river$FD$toRN[j]==0){
      distanz <- distanz + river$FD$leng[j]
      j <- river$FD$downNode[j]}
    distanz <- distanz + river$FD$leng[j]
    RNnode <- river$FD$toRN[j]
  }

  AGnode <- river$RN$toAGReach[RNnode]
  Xnew <- river$RN$X[RNnode]; Ynew <- river$RN$Y[RNnode]

  if (showPlot){
    Xmin <- min(X,Xnew)-20*river$cellsize; Xmax <- max(X,Xnew)+20*river$cellsize;
    Ymin <- min(Y,Ynew)-20*river$cellsize; Ymax <- max(Y,Ynew)+20*river$cellsize

    if (is.null(xlim)){xlim <- c(Xmin, Xmax)}
    if (is.null(ylim)){ylim <- c(Ymin, Ymax)}

    theme <- numeric(river$AG$nNodes); theme[AGnode] <- 1
    plot(river,theme, discreteLevels = T,
         colPalette = colorRampPalette(c("blue2","orange")),
         xlim=xlim, ylim=ylim)

    xy_lim <- par("usr")

    legend(x=xy_lim[1]+river$cellsize, y=xy_lim[4], legend=c("Original site","RN node","AG node"),
           col=c("red","black","orange"), pch=c(15,20,NA),lty=c(0,0,1))
    points(X,Y,pch=15,col="red")
    points(Xnew,Ynew,pch=20,col="black")

  }

  explist <- vector("list")
  explist[["FDnode"]] <- indFD
  explist[["distance"]] <- distanz
  explist[["AGnode"]] <- AGnode
  explist[["RNnode"]] <- RNnode

  invisible(explist)
}
