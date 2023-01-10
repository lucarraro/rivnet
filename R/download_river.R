download_river <- function(outlet=c(735010, 261530),
                           ext=c(700000, 770000, 220000, 270000),
                           EPSG=21781,
                           zoom=9){

test_dir <- withr::local_tempdir() # temporary directory storing intermediary files created by TauDEM

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

out_shp <-  shp.point(x_outlet,y_outlet,file.path(test_dir,"ApproxOutlet")) # create shp for approximate outlet

# Move outlet to stream
out_moved.shp <- traudem::taudem_moveoutletstostream(out_p, out_src, outlet_file = out_shp,
                                                     output_moved_outlets_file = file.path(test_dir,"Outlet.shp"), quiet = TRUE)

# Contributing area upstream of outlet
out_ssa <- traudem::taudem_aread8(out_p, outlet_file = out_moved.shp, quiet = TRUE)

# Derive spatRaster from TauDEM output for subsequent elaboration
fel <- rast(out_fel) # pit-filled DEM
p <- rast(out_p) # DB flow direction map
ssa <- rast(out_ssa) # contributing area within the catchment contour

out <- list(fel=fel, p=p, ssa=ssa)
invisible(out)
}

shp.point <- function(x,y,sname="shape"){ # function to create point shapefile given coordinates
  n <- length(x)
  dd <- data.frame(Id=1:n,X=x,Y=y)
  ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
  invisible(paste0(sname,".shp"))
}
