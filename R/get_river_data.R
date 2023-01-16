get_river_data <- function(outlet=c(735010, 261530),
                           ext=c(700000, 770000, 220000, 270000),
                           EPSG=21781,
                           zoom=9,
                           DEM=NULL,
                           as.channelNetwork=TRUE,
                           as.rast=FALSE,
                           filename=NULL){

  if (as.rast==T & is.null(filename)){
    stop("filename cannot be NULL if as.rast = TRUE.")
  }

  test_dir <- withr::local_tempdir() # temporary directory storing intermediary files created by TauDEM

  if (is.null(DEM)){
  # Input data for the river Thur ####
  #EPSG <- 21781	# EPSG code corresponding to CH1903 projection. See also http://www.epsg-registry.org
  x_outlet <- outlet[1] #735010 # outlet x coordinate
  y_outlet <- outlet[2] #261530 # outlet y coordinate
  x_ll <- ext[1] #700000 # x coordinate of the lower-left (SW) limit of the region
  x_tr <- ext[2] #770000 # x coordinate of the top-right (NE) limit of the region
  y_ll <- ext[3] #220000 # y coordinate of the lower-left (SW) limit of the region
  y_tr <- ext[4] #270000 # y coordinate of the top-right (NE) limit of the region
  # zoom <- 9 # corresponds to cellsize ~= 100 m

  # use elevatr to download DEM
  loc.df <- data.frame(x=c(x_ll, x_tr), y=c(y_ll,y_tr))
  r <- rast()
  crs(r) <- paste0("epsg:",EPSG)
  crs_str <- crs(r)

  z <- get_elev_raster(locations = loc.df, prj = crs_str, z=zoom, verbose=F) # call elevatr
  z <- rast(z) # transform into spatRaster object
  z <- crop(z, ext(x_ll, x_tr, y_ll, y_tr)) # crop spatRaster to the region of interest (to save RAM)
  z <- classify(z, matrix(c(NA,NA,0),1,3)) # all pixels with elev=NA are set to 0. Then the pit remove algorithm will take care of them
  }
  else {z <- rast(DEM)}

  writeRaster(z,filename=file.path(test_dir, "DEM.tif")) # write elevation raster to temporary directory

  # apply TauDEM functions
  # Remove pits
  out_fel <- traudem::taudem_pitremove(file.path(test_dir, "DEM.tif"), quiet = TRUE)

  # D8 flow directions
  out_p <- traudem::taudem_d8flowdir(out_fel, quiet = TRUE)
  out_p <- out_p$output_d8flowdir_grid # file path of flow direction file

  # Contributing area
  out_ad8 <- traudem::taudem_aread8(out_p, quiet = TRUE)

  # Threshold
  out_src <- traudem::taudem_threshold(out_ad8, quiet = TRUE)

  p.sf <- sf::st_as_sf(data.frame(x = x_outlet,y= y_outlet), coords = c("x", "y"), crs = EPSG)
  out_shp <- file.path(test_dir,"ApproxOutlet.shp")
  sf::st_write(p.sf, out_shp, driver="ESRI Shapefile", quiet=T)

  #out_shp <-  shp.point(x_outlet,y_outlet,file.path(test_dir,"ApproxOutlet")) # create shp for approximate outlet

  # Move outlet to stream
  out_moved.shp <- traudem::taudem_moveoutletstostream(out_p, out_src, outlet_file = out_shp,
                                                       output_moved_outlets_file = file.path(test_dir,"Outlet.shp"), quiet = TRUE)

  # Contributing area upstream of outlet
  out_ssa <- traudem::taudem_aread8(out_p, outlet_file = out_moved.shp, quiet = TRUE)

  # Derive spatRaster from TauDEM output for subsequent elaboration
  fel <- rast(out_fel) # pit-filled DEM
  p <- rast(out_p) # DB flow direction map
  ssa <- rast(out_ssa) # contributing area within the catchment contour

  if (as.rast==T){
    rr <- c(fel, p, ssa) # export this if you can yield the source data
    names(rr) <- c("pit-filledDEM","flowDir","drainageArea")
    writeRaster(rr,filename)
  }

  if (as.channelNetwork==TRUE){
    ncols <- ncol(ssa)
    nrows <- nrow(ssa)
    ssa_val <- values(ssa)
    flowDir <- values(p)

    Nnodes_FD <- sum(!is.na(ssa_val))
    FD_to_DEM <- which(!is.na(ssa_val))

    X_FD <- xFromCell(ssa,1:ncell(ssa))[FD_to_DEM]
    Y_FD <- yFromCell(ssa,1:ncell(ssa))[FD_to_DEM]
    Z_FD <- values(fel)[FD_to_DEM]
    cellsize <- sqrt(prod(res(fel)))
    A_FD <- cellsize^2*ssa_val[FD_to_DEM]

    Length <- cellsize*(1 + as.numeric(flowDir %% 2 ==0)*(sqrt(2)-1))
    Length_FD <- Length[FD_to_DEM]

    nNodes_FD <- length(FD_to_DEM)
    W_FD <- spam(0,nNodes_FD,nNodes_FD)
    ind <- matrix(0,nNodes_FD,2)
    downNode_FD <- numeric(nNodes_FD)
    Slope_FD <-  numeric(nNodes_FD)

    k <- 1
    for (i in 1:nNodes_FD){
      mov <- neigh(flowDir[FD_to_DEM[i]])
      d <- which(FD_to_DEM==(FD_to_DEM[i]+mov[1]+mov[2]*ncols)) # indices from top-left corner to the right, then next row...
      if (length(d)!=0){
        ind[k, ] <- c(i,d)
        k <- k + 1
        Slope_FD[i] <- (Z_FD[i]-Z_FD[d])/Length_FD[i]
      }
    }
    ind <- ind[-k, ]
    downNode_FD[ind[,1]] <- ind[,2]
    W_FD[ind] <- 1
    rm(ind)
    Outlet_FD <- which(downNode_FD==0)

    newW <- new("spam")
    slot(newW, "entries", check = FALSE) <- W_FD@entries[-1]
    slot(newW, "colindices", check = FALSE) <- W_FD@colindices[-1]
    slot(newW, "rowpointers", check = FALSE) <- c(W_FD@rowpointers[1],W_FD@rowpointers[-1]-W_FD@rowpointers[1])
    slot(newW, "dimension", check = FALSE) <- W_FD@dimension
    W_FD <- newW

    pl <- initial_permutation(downNode_FD)
    pl <- as.integer(pl$perm)

    ssa_cont <- classify(ssa,matrix(c(NA,NA,-1e6,1,Inf,1e6),2,3,byrow=T))
    cont <- as.contour(ssa_cont,levels=c(0,1e6))
    cont <- subset(cont, cont$level==0) # pick most external contour
    XContour <-  crds(cont)[,1]
    YContour <-  crds(cont)[,2]

    CM <- list(A=max(A_FD,na.rm=T), XContour=list(list(XContour)), YContour=list(list(YContour)),
               XContourDraw=list(list(XContour)), YContourDraw=list(list(YContour)))
    FD <- list(A=A_FD,W=W_FD,downNode=downNode_FD,X=X_FD,Y=Y_FD,nNodes=Nnodes_FD,outlet=Outlet_FD,perm=pl,
               Z=Z_FD,slope=Slope_FD,leng=Length_FD,XDraw=X_FD,YDraw=Y_FD,toCM=1+numeric(Nnodes_FD))

    river <- list(FD=FD,CM=CM,dimX=ncols,dimY=nrows,cellsize=cellsize,nOutlet=1,periodicBoundaries=FALSE,
                  expEnergy=NULL,coolingRate=NULL,typeInitialState=NULL,nIter=NULL,initialNoCoolingPhase=NULL,
                  energyInit=NULL)

    invisible(river)
  }
}

# shp.point <- function(x,y,sname="shape"){ # function to create point shapefile given coordinates
#   n <- length(x)
#   dd <- data.frame(Id=1:n,X=x,Y=y)
#   ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
#   ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
#   write.shapefile(ddShapefile, sname, arcgis=T)
#   invisible(paste0(sname,".shp"))
# }

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



