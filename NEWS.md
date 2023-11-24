# rivnet 0.3.3.9000

## Minor changes

- `extract_river`: a warning is thrown whenever `EPSG` identifies a geographic projection system.
- `extract_river`: option `showPlot = TRUE` shows an additional plot with a zoom-in in the outlet 
area, which can help diagnose issues in catchment delineation.
- `points_colorscale`: Argument `...` added.

# rivnet 0.3.3

## Minor changes

- `DESCRIPTION` updated with DOI of published article.
- `CITATION` added. 

# rivnet 0.3.2

## Minor changes

- `extract_river`: option `args_get_elev_raster` is added.

# rivnet 0.3.1

## Major changes

- `points_colorscale`: function added.

# rivnet 0.3.0

## Major changes

- `river_to_AEM`: function added.

# rivnet 0.2.1

## Minor changes

- `extract_river`: `curl::has_internet()` check by `elevatr::get_elev_raster` is overridden.
- `extract_river`: crs attributed when outlet shapefile is written in temporary directory.
- `locate_site`: `par` is exported when `showPlot = TRUE`. Example updated.

# rivnet 0.2.0

## Major changes

- `"[["`, `"[[<-"` methods for `river` class defined and examples added. 
- `path_velocities_river`: C++ implementation added.
- Dependency on `Rcpp` added.
- `covariate_river`: argument `overwrite` added.
- `locate_site`: site coordinates can be passed as a list. Example added.

## Minor changes

- Info on installation from CRAN updated in vignette.
- `inst/extdata/temperature.tif` changed.
- `extract_river`: a warning is thrown when the WGS84 (`EPSG = 4326`) projection system is selected.

## Bugs solved

- `extract_river`: multiple catchment shapes are shown when `showPlot = TRUE` 
and multiple, non-neighboring catchments are extracted.
- `covariate_river`: fixed bug for raster file containing `NA` values.
