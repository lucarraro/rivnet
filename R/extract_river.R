extract_river <- function(outlet,
                          EPSG=NULL,
                          ext=NULL,
                          z=NULL,
                          DEM=NULL,
                          as.river=TRUE,
                          as.rast=FALSE,
                          filename=NULL,
                          showPlot=FALSE,
                          threshold_parameter=1000,
                          n_processes=1,
                          displayUpdates=0,
                          src="aws",
                          args_get_elev_raster=list()){

  if (!is.null(EPSG)){
    if ( is.null(st_crs(EPSG)$units)){
      warning("You are using a geographic coordinate system.
    It is recommended that you use a projected coordinate system.")}}

  if (!is.null(DEM)){
    elev <- rast(DEM)
    ext <- as.vector(ext(elev))}

  if (as.rast==T & is.null(filename)){
    stop("filename cannot be NULL if as.rast = TRUE.")
  }

  if (ext[2] < ext[1] | ext[4] < ext[3]){
    stop('ext has wrong format. It should be provided as c(xmin, xmax, ymin, ymax).')
  }

  if (is.vector(outlet)) {outlet <- matrix(outlet, nrow=length(outlet)/2, ncol=2, byrow=T)}

  if (any(outlet[,1]>ext[2] | outlet[,1]<ext[1] | outlet[,2]>ext[4] | outlet[,2]<ext[3]) ){
    stop('outlet coordinates are beyond the region defined by ext')
  }

  if ("x" %in% names(outlet)){x_outlet <- outlet$x} else {x_outlet <- outlet[,1]}
  if ("y" %in% names(outlet)){y_outlet <- outlet$y} else {y_outlet <- outlet[,2]}

  #test_dir <- withr::local_tempdir() # temporary directory storing intermediary files created by TauDEM
  test_dir <- tempdir()

  if (displayUpdates==2) {quiet=FALSE} else {quiet=TRUE}

  t0 <- Sys.time()
  if (is.null(DEM)){

    # use elevatr to download DEM
    loc.df <- data.frame(x=c(ext[1], ext[2]), y=c(ext[3],ext[4]))
    r <- rast()
    crs(r) <- paste0("epsg:",EPSG)
    crs_str <- crs(r)

    # override curl::has_internet() check by get_elev_raster
    op <- get("has_internet_via_proxy", environment(curl::has_internet)) # old value
    # check for internet
    np <- !is.null(curl::nslookup("r-project.org", error = FALSE))
    assign("has_internet_via_proxy", np, environment(curl::has_internet))
    on.exit(assign("has_internet_via_proxy", op, environment(curl::has_internet)))

    if (is.null(args_get_elev_raster$locations)) args_get_elev_raster$locations=loc.df
    if (is.null(args_get_elev_raster$prj)) args_get_elev_raster$prj=crs_str
    if (is.null(args_get_elev_raster$z)) args_get_elev_raster$z=z
    if (is.null(args_get_elev_raster$verbose)) args_get_elev_raster$verbose=(displayUpdates>0)
    if (is.null(args_get_elev_raster$clip)) args_get_elev_raster$clip="bbox"
    if (is.null(args_get_elev_raster$src)) args_get_elev_raster$src=src

    if (quiet){
      elev <- suppressMessages(do.call(get_elev_raster, args_get_elev_raster))
    } else {
      elev <- do.call(get_elev_raster, args_get_elev_raster)}
    elev <- rast(elev) # transform into spatRaster object
    elev <- classify(elev, matrix(c(NA,NA,0),1,3)) # all pixels with elev=NA are set to 0. Then the pit remove algorithm will take care of them
  }

  writeRaster(elev, filename=file.path(test_dir, "DEM.tif"), overwrite=TRUE) # write elevation raster to temporary directory
  t1 <- Sys.time()

  # apply TauDEM functions
  # Remove pits
  if (displayUpdates==1) message('Remove pits...\n',appendLF=FALSE)
  out_fel <- traudem::taudem_pitremove(file.path(test_dir, "DEM.tif"),
                                       n_processes = n_processes,
                                       quiet = quiet)

  # D8 flow directions
  if (displayUpdates==1) message('D8 flow directions...\n',appendLF=FALSE)
  out_p <- traudem::taudem_d8flowdir(out_fel,
                                     n_processes = n_processes,
                                     quiet = quiet)
  out_p <- out_p$output_d8flowdir_grid # file path of flow direction file

  # Contributing area
  if (displayUpdates==1) message('Contributing areas...\n',appendLF=FALSE)
  out_ad8 <- traudem::taudem_aread8(out_p,
                                    n_processes = n_processes,
                                    quiet = quiet)

  # Threshold
  if (displayUpdates==1) message('Stream definition by threshold...\n',appendLF=FALSE)
  out_src <- traudem::taudem_threshold(out_ad8,
                                       threshold_parameter = threshold_parameter,
                                       n_processes = n_processes,
                                       quiet = quiet)

  if (!is.null(EPSG)){
    p.sf <- sf::st_as_sf(data.frame(x = x_outlet,y= y_outlet), coords = c("x", "y"), crs=EPSG) # crs=EPSG
  } else {
    p.sf <- sf::st_as_sf(data.frame(x = x_outlet,y= y_outlet), coords = c("x", "y"))
  }

  out_shp <- file.path(test_dir,"ApproxOutlet.shp")
  sf::st_write(p.sf, out_shp, driver="ESRI Shapefile", quiet=quiet, append=FALSE)


  # Move outlet to stream
  if (displayUpdates==1) message('Move outlet to stream...\n',appendLF=FALSE)
  out_moved.shp <- traudem::taudem_moveoutletstostream(out_p, out_src, outlet_file = out_shp,
                                                       output_moved_outlets_file = file.path(test_dir,"Outlet.shp"),
                                                       n_processes = n_processes,
                                                       quiet = quiet)


  # Contributing area upstream of outlet
  if (displayUpdates==1) message('Contributing area upstream of outlet...\n',appendLF=FALSE)
  out_ssa <- traudem::taudem_aread8(out_p, output_contributing_area_grid = file.path(test_dir,"ssa.tif"),
                                    n_processes = n_processes,
                                    outlet_file = out_moved.shp, quiet = quiet)

  # Derive spatRaster from TauDEM output for subsequent elaboration
  fel <- rast(out_fel) # pit-filled DEM
  p <- rast(out_p) # DB flow direction map
  ad8 <- rast(out_ad8) # contributing area for the whole region
  ssa <- rast(out_ssa) # contributing area within the catchment contour

  shp <- sf::st_read(out_moved.shp,quiet=T)
  out_moved <- sf::st_coordinates(shp)
  cellsize <- sqrt(prod(res(fel)))
  no.cells <- as.numeric(as.matrix(extract(ad8, out_moved)))
  if (length(no.cells)==1){
    outlet <- matrix(outlet,nrow=1,ncol=2,byrow=T)
    out_moved <- matrix(out_moved,nrow=1,ncol=2,byrow=T)
  }

  if (showPlot==T){
    oldpar <- par(no.readonly = TRUE)
    if (length(no.cells)>1){
      par(mfrow=c(1,2))
      on.exit(par(oldpar))
    }
    for (ind_out in 1:length(no.cells)){
    plot(ad8,col=hcl.colors(1000),
         xlim=out_moved[ind_out, 1]+cellsize*c(-50,50),
         ylim=out_moved[ind_out, 2]+cellsize*c(-50,50))
    points(outlet[ind_out, 1],outlet[ind_out, 2], pch=22,bg="red")
    points(out_moved[ind_out, 1],out_moved[ind_out, 2],
           pch=21 ,bg="magenta",cex=0.8)
    legend(out_moved[ind_out, 1]-50*cellsize, out_moved[ind_out, 2]+50*cellsize,
           pt.bg=c("red","magenta"), pch=c(22,21), pt.cex=c(1,0.8),cex=c(0.75,0.75),
           legend=c("Original","Moved"), title="Outlet")
    title(sprintf("Catchment contains %d pixels",no.cells[ind_out]))
    }
    par(oldpar)
  }

  t2 <- Sys.time()
  # create catchment contour
  ssa_cont <- classify(ssa,matrix(c(NA,NA,-1e6,1,Inf,1e6),2,3,byrow=T))
  cont <- as.contour(ssa_cont,levels=c(0,1e6))
  cont <- subset(cont, cont$level==0) # pick most external contour
  XC <-  crds(cont)[,1]
  YC <-  crds(cont)[,2]
  # patch for irregular contour: cut part of a contour that is unconnected to the rest
  ind_cut <- which(abs(diff(XC))>10*cellsize | abs(diff(YC))>10*cellsize)
  XContour <- YContour <- list()
  ind_cut <- c(0,ind_cut,length(XC))
  for (i in 1:(length(ind_cut)-1)){
    XContour[[i]] <- XC[(ind_cut[i]+1):ind_cut[i+1]]
    YContour[[i]] <- YC[(ind_cut[i]+1):ind_cut[i+1]]
  }

  if (showPlot==T){
    plot(ad8,col=hcl.colors(1000))
    for (i in 1:length(XContour)){
      lines(XContour[[i]],YContour[[i]],col="magenta")}
    title(sprintf('Max drainage area: %.2e m2',max(values(ssa)*prod(res(ssa)),na.rm=T)))
  }

  if (as.rast==T){
    rr <- c(fel, p, ad8, ssa) # export this if you can yield the source data
    names(rr) <- c("fel","p","ad8", "ssa")
    writeRaster(rr,filename,overwrite=T)
  }

  if (as.river==TRUE){
    if (displayUpdates>0){message("Creation of river object...      \r", appendLF = FALSE)}
    ncols <- ncol(ssa)
    nrows <- nrow(ssa)
    ssa_val <- values(ssa)
    flowDir <- values(p)

    Nnodes_FD <- sum(!is.na(ssa_val))
    FD_to_DEM <- which(!is.na(ssa_val))
    DEM_to_FD <- numeric(max(FD_to_DEM))
    DEM_to_FD[FD_to_DEM] <- 1:Nnodes_FD

    X_FD <- xFromCell(ssa,1:ncell(ssa))[FD_to_DEM]
    Y_FD <- yFromCell(ssa,1:ncell(ssa))[FD_to_DEM]
    Z_FD <- values(fel)[FD_to_DEM]
    A_FD <- cellsize^2*ssa_val[FD_to_DEM]
    xllcorner <- min(xFromCell(ssa,1:ncell(ssa)))
    yllcorner <- min(yFromCell(ssa,1:ncell(ssa)))

    Length <- cellsize*(1 + as.numeric(flowDir %% 2 ==0)*(sqrt(2)-1))
    Length_FD <- Length[FD_to_DEM]

    nNodes_FD <- length(FD_to_DEM)
    W_FD <- spam(0,nNodes_FD,nNodes_FD)
    ind <- matrix(0,nNodes_FD,2)
    downNode_FD <- numeric(nNodes_FD)
    Slope_FD <-  numeric(nNodes_FD)

    downNode_rev <- vector(nNodes_FD,mode="list")

    k <- 1
    for (i in 1:nNodes_FD){
      mov <- neigh(flowDir[FD_to_DEM[i]])
      d <- DEM_to_FD[FD_to_DEM[i]+mov[1]+mov[2]*ncols]# indices from top-left corner to the right, then next row...
      #d <- which(FD_to_DEM==(FD_to_DEM[i]+mov[1]+mov[2]*ncols)) # slow alternative
      if (is.na(d)){d=0}
      if (d!=0){
        ind[k, ] <- c(i,d)
        k <- k + 1
        Slope_FD[i] <- (Z_FD[i]-Z_FD[d])/Length_FD[i]
        downNode_rev[[d]] <- c(downNode_rev[[d]],i)
      }
      # if (!quiet){
      #   if ((i %% round(nNodes_FD*0.01))==0){
      #   message(sprintf("Creation of river object... %.1f%% \r",i/(1.001*nNodes_FD)*100), appendLF = FALSE)}}
    }
    ind <- ind[-(k:nNodes_FD), ]
    downNode_FD[ind[,1]] <- ind[,2]
    W_FD[ind] <- 1
    rm(ind)
    Outlet_FD <- which(downNode_FD==0)
    if (length(Outlet_FD) != dim(outlet)[1]){
      stop("The number of identified outlets is not equal to the number of inputted outlets.
       Some of the extracted catchments might be nested.
       If this is the case, run extract_river separately for each catchment.")
    }
    # reassign outlet order
    if (length(Outlet_FD)>1){
      ind_out <- numeric(length(Outlet_FD))
      for (i in 1:length(Outlet_FD)){
        dist <- sqrt((X_FD[Outlet_FD[i]]-out_moved[,1])^2+(Y_FD[Outlet_FD[i]]-out_moved[,2])^2)
        ind_out[i] <- which(dist==min(dist))
      }
      Outlet_FD <- Outlet_FD[ind_out]
    }

    newW <- new("spam")
    slot(newW, "entries", check = FALSE) <- W_FD@entries[-1]
    slot(newW, "colindices", check = FALSE) <- W_FD@colindices[-1]
    slot(newW, "rowpointers", check = FALSE) <- c(W_FD@rowpointers[1],W_FD@rowpointers[-1]-W_FD@rowpointers[1])
    slot(newW, "dimension", check = FALSE) <- W_FD@dimension
    W_FD <- newW

    #pl <- initial_permutation_rev(downNode_rev, Outlet_FD)
    pl <- init_perm_rev_cpp(downNode_rev, Outlet_FD)
    pl <- as.integer(pl$perm)

    toCM <- numeric(Nnodes_FD)
    for (i in 1:length(Outlet_FD)){
      tmp <- which(pl==Outlet_FD[i])
      toCM[pl[(tmp - A_FD[Outlet_FD[i]]/cellsize^2 + 1):tmp]] <- i
    }

    # draw contours of single catchments if length(Outlet_FD)>1
    if (length(Outlet_FD)>1){
      XContour <- YContour <- list()
      for (i in 1:length(Outlet_FD)){
        ssa_tmp <- ssa
        mask <- NA*numeric(length(values(ssa_tmp)))
        mask[FD_to_DEM[which(toCM==i)]] <- 1
        values(ssa_tmp) <- values(ssa_tmp)*mask

        ssa_cont <- classify(ssa_tmp,matrix(c(NA,NA,-1e6,1,Inf,1e6),2,3,byrow=T))
        cont <- as.contour(ssa_cont,levels=c(0,1e6))
        cont <- subset(cont, cont$level==0) # pick most external contour
        XCont <-  crds(cont)[,1]
        YCont <-  crds(cont)[,2]
        # patch for irregular contour: cut part of a contour that is unconnected to the rest
        ind_cut <- which(abs(diff(XCont))>10*cellsize | abs(diff(YCont))>10*cellsize)
        if (length(ind_cut)>0){
          XCont <- XCont[1:ind_cut]
          YCont <- YCont[1:ind_cut]
        }
        XContour <- c(XContour, list(list(XCont)))
        YContour <- c(YContour, list(list(YCont)))
      }
    } else {
      XContour <- list(XContour)
      YContour <- list(YContour)
    }

    if (displayUpdates>0){message("Creation of river object... 100.0% \n", appendLF = FALSE)}

    CM <- list(A=A_FD[Outlet_FD], XContour=XContour, YContour=YContour,
               XContourDraw=XContour, YContourDraw=YContour)
    FD <- list(A=A_FD,W=W_FD,downNode=downNode_FD,X=X_FD,Y=Y_FD,nNodes=Nnodes_FD,outlet=Outlet_FD,perm=pl,
               Z=Z_FD,slope=Slope_FD,leng=Length_FD,XDraw=X_FD,YDraw=Y_FD,toCM=toCM, toDEM=FD_to_DEM)

    river <- list(FD=FD,CM=CM,dimX=ncols,dimY=nrows,cellsize=cellsize,nOutlet=length(Outlet_FD),periodicBoundaries=FALSE,
                  xllcorner=xllcorner, yllcorner=yllcorner)

    river_S4 <- new("river")
    fieldnames <- names(river)
    for (i in 1:length(fieldnames)){slot(river_S4, fieldnames[i]) <- river[[fieldnames[i]]]}
  }
  t3 <- Sys.time()
  if (displayUpdates>0){
    message("extract_river has finished. \n",appendLF = FALSE)
    message(sprintf("Time for DEM download: %.1f s \n",difftime(t1,t0,units="secs")),appendLF = FALSE)
    message(sprintf("Time for TauDEM processing: %.1f s \n",difftime(t2,t1,units="secs")),appendLF = FALSE)
    message(sprintf("Time for creation of river object: %.1f s \n",difftime(t3,t2,units="secs")),appendLF = FALSE)
  }

  unlink(file.path(test_dir,"*"))

  if (as.river==TRUE){invisible(river_S4)}
}

setClass("river",
         slots= c(FD="list", dimX="numeric", dimY="numeric", cellsize="numeric", nOutlet="numeric",
                  periodicBoundaries="logical", expEnergy="numeric", coolingRate="numeric",
                  typeInitialState="character", nIter="numeric", initialNoCoolingPhase="numeric",
                  energy="numeric", exitFlag="numeric", N4="list",N8="list", nIterSequence="numeric",
                  energyInit="numeric", xllcorner="numeric", yllcorner="numeric",
                  CM="list",RN="list",AG="list",OptList="list",SC="list",thrA="numeric",
                  slope0="numeric", zMin="numeric",streamOrderType="character",maxReachLength="numeric",
                  widthMax="numeric",depthMax="numeric",velocityMax="numeric",expWidth="numeric",expDepth="numeric",expVelocity="numeric"))

setMethod("$","river",
          function(x,name){
            slot(x,name)
          })


setMethod("show", "river",
          function(object){
            isOCN <- length(object$coolingRate) > 0
            isElev <- length(object$CM) > 0
            isAggr <- length(object$RN) > 0
            isPath <- length(object$AG$downstreamPath) > 0
            isRivG <- length(object$AG$width) > 0

            cat("Class         : river \n")
            if (isOCN){
              if (object$typeInitialState=="custom"){
                cat("Type          : Optimal Channel Network (general contour) \n")
              } else {
                cat("Type          : Optimal Channel Network \n")}
            } else {
              cat("Type          : Real river \n")
            }
            cat(sprintf("No. FD nodes  : %d \n", object$FD$nNodes))
            cat(sprintf('Dimensions    : %d x %d \n',object$dimX,object$dimY))
            cat(sprintf('Cell size     : %.2f \n',object$cellsize))
            cat("Has elevation :",isElev,"\n" )
            if (isElev){
              cat("Aggregated    :",isAggr,"\n" )
              if (isAggr){
                cat(sprintf("   Threshold area  : %.2f \n",object$thrA))
                cat(sprintf("   Max reach length: %.2f \n",object$maxReachLength))
                cat(sprintf("   No. RN nodes    : %d \n",object$RN$nNodes))
                cat(sprintf("   No. AG nodes    : %d \n",object$AG$nNodes))
                cat(        "   Has paths       : ",isPath,"\n")
                cat(        "   River geometry  : ",isRivG,"\n")
              }}
          })


setMethod("names",signature=c(x="river"),
          function(x){
            slotNames(x)
          })

setMethod("$<-",signature=c(x="river"),
          function(x,name,value){
            slot(x, name) <- value
            return(x)
          })

setMethod("[[",c("river","character","missing"),
          function(x,i){
            slot(x,i)
          })

setMethod("[[<-",signature=c("river","character","missing"),
          function(x,i,value){
            slot(x, i) <- value
            return(x)
          })

setMethod("plot", signature(x="river",y="missing"),
          function(x, type, ...){
            if (missing(type)) {type <- "RN"}
            if (type=="elev2D") {OCNet::draw_elev2D_OCN(x, ...)}
            else if (length(x$AG) > 0) {
              if (type=="SC" | type=="subcatchments"){
                OCNet::draw_subcatchments_OCN(x, ...)
              } else {OCNet::draw_thematic_OCN(x, ...)}
            } else if (length(x$CM) > 0) {
              OCNet::draw_contour_OCN(x, ...)
            } else {
              OCNet::draw_simple_OCN(x, ...)}})

setMethod("plot", signature(x="numeric",y="river"),
          function(x, y, type, ...){
            if (missing(type)) type <- "RN"
            if (isTRUE(type=="SC" | type=="subcatchments")) {
              OCNet::draw_subcatchments_OCN(y, x, ...)
            } else { OCNet::draw_thematic_OCN(y, x, ...)}})

setMethod("plot", signature(x="river",y="numeric"),
          function(x, y, type, ...){
            if (missing(type)) type <- "RN"
            if (isTRUE(type=="SC" | type=="subcatchments")) {
              OCNet::draw_subcatchments_OCN(x, y, ...)
            } else { OCNet::draw_thematic_OCN(y, x, ...)}
          })


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

# deprecated
initial_permutation_rev <- function(downNode_rev, Outlet){
  nNodes <- length(downNode_rev)
  NodesToExplore <- Outlet # start from outlets
  reverse_perm <- numeric(nNodes) # build permutation vector from outlets to headwaters, then flip it
  k <- 0
  while (length(NodesToExplore)>0){ # continue until all the network has been explored
    k <- k + 1
    node <- NodesToExplore[1] # explore a node
    reverse_perm[k] <- node # assign position in the permutation vector
    NodesToExplore <- NodesToExplore[-1] # remove explored node
    UpNodes <- downNode_rev[[node]] # find nodes upstream of node
    while (length(UpNodes)>0){ # continue upstream until a headwater is found
      k <- k + 1
      node <- UpNodes[1] # explore first upstream node
      reverse_perm[k] <- node
      if (length(UpNodes)>1){ # if there is a bifurcation upstream, add the other upstream connections at the top of NodesToExplore
        NodesToExplore <- c(UpNodes[2:length(UpNodes)],NodesToExplore)
      }
      UpNodes <- downNode_rev[[node]]
    }
  }
  perm <- reverse_perm[nNodes:1] # flip permutation
  OutList = list(perm=perm,noDAG=0)
  invisible(OutList)
}

