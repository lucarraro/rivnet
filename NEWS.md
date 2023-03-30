# rivnet 0.1.0.9000

## Major changes

- `"[["`, `"[[<-"` methods for `river` class defined and examples added. 
- `path_velocities_river`: C++ implementation added.
- Dependency on `Rcpp` added.
- `covariate_river`: argument `overwrite` added.

## Minor changes

- Info on installation from CRAN updated in vignette.
- `inst/extdata/temperature.tif` changed.
- `extract_river`: a warning is thrown when the WGS84 (`EPSG = 4326`) projection system is selected.

## Bugs solved

- `extract_river`: multiple catchment shapes are shown when `showPlot = TRUE` 
and multiple, non-neighboring catchments are extracted.