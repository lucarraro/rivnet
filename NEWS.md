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
