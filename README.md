<img align="right" width="250" src="man/figures/rivnet_logo.png">

# rivnet

An R-package allowing seamless extraction of river networks from Digital Elevation Models data.

[![CRAN](http://www.r-pkg.org/badges/version/rivnet)](http://CRAN.R-project.org/package=rivnet)

## What this package does 

- Analysis of Digital Elevation Models, either provided externally or downloaded from open source repositories (in the latter case, interfacing with the `elevatr` package). 
- Extraction of river networks is performed via TauDEM's D8 flow direction algorithm (interfacing with the `traudem` package). 
- Resulting river networks are compatible with functions from the `OCNet` package, and can be plotted and analyzed accordingly (and are thus also compatible with packages `igraph` and `SSN`). 
- Distances, areas, subcatchments, slopes and so on can be computed. The obtained river networks can be used for further hydrological and ecological modelling studies.

## A minimal working example

Extract the river Wigger (Switzerland) from an externally provided DEM raster file. Outlet coordinates are in the CH1903/LV03 projected coordinate system (i.e., the same as the DEM file):

```
 fp <- system.file("extdata/wigger.tif", package = "rivnet")
 r <- extract_river(outlet = c(637478, 237413),
	                  DEM = fp)
````

The same result can be obtained by downloading DEM data via `elevatr`:

```
r <- extract_river(outlet = c(637478, 237413),
	                  EPSG = 21781, #CH1903/LV03 coordinate system
	                  ext = c(6.2e5, 6.6e5, 2e5, 2.5e5),
	                  z = 8)

````

## Installation issues

`rivnet` depends on `traudem`, whose installation might require some caution depending on your operating system. Please read [its documentation](https://lucarraro.github.io/traudem/) carefully.

## Author

Luca Carraro (maintainer)
