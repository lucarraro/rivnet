locate_site <- function(X,Y=NULL,river,
                        euclidean=TRUE,
                        showPlot=FALSE,
                        xlim=NULL,
                        ylim=NULL){

  if (missing(river)){
      if (is(Y,"river")){river <- Y
      } else stop("river object must be provided.")}

  if (length(river$RN$X)==0){
    stop('Missing fields in river. You should run aggregate_river prior to locate_site.')
  }

  if (is.list(X)){
    if ("x" %in% names(X)){
      x <- X$x
    } else if ("X" %in% names(X)){
      x <- X$X
    } else stop("If X is a data frame, it must contain object x or X")
    if ("y" %in% names(X)){
      y <- X$y
    } else if ("Y" %in% names(X)){
      y <- X$Y
    } else stop("If X is a data frame, it must contain object y or Y")
    X <- x; Y <- y
  } else if(is.numeric(X) & is.null(Y)){
    stop("If X is numeric, Y cannot be NULL.")
  }

  if (length(X)>1 | length(Y)>1) stop("Coordinates cannot have length > 1.")

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
    #oldpar <- par(no.readonly = TRUE)
    #on.exit(par(oldpar))

    Xmin <- min(X,Xnew)-20*river$cellsize; Xmax <- max(X,Xnew)+20*river$cellsize;
    Ymin <- min(Y,Ynew)-20*river$cellsize; Ymax <- max(Y,Ynew)+20*river$cellsize

    if (is.null(xlim)){xlim <- c(Xmin, Xmax)}
    if (is.null(ylim)){ylim <- c(Ymin, Ymax)}

    theme <- numeric(river$AG$nNodes); theme[AGnode] <- 1
    plot(river,theme, discreteLevels = T,
         colPalette = colorRampPalette(c("blue2","orange")),
         xlim=xlim, ylim=ylim)
    xy_lim <- par("usr")

    legend(x=xy_lim[1]+river$cellsize, y=xy_lim[4], legend=c("Original site","RN node","AG reach"),
           col=c("red","black","orange"), pch=c(15,20,NA),lty=c(0,0,1))
    points(X,Y,pch=15,col="red")
    points(Xnew,Ynew,pch=20,col="black")

    parfig <- par() # useless?
  }

  explist <- vector("list")
  explist[["FDnode"]] <- indFD
  explist[["distance"]] <- distanz
  explist[["AGnode"]] <- AGnode
  explist[["RNnode"]] <- RNnode
  if (showPlot) explist[["par"]] <- parfig

  invisible(explist)
}
