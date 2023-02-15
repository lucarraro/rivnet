\name{covariate_river}
\alias{covariate_river}
\title{Attribute covariates to nodes of a river network}
\usage{
covariate_river(x, river, categorical = TRUE)
}
\description{
Attributes covariate values from raster files to subcatchments of a river network. Both local and upstream-averaged 
covariate values are calculated.
}

\arguments{
\item{x}{SpatRaster object (obtained via \code{rast}) containing categorical or continuous variables from which covariates
are computed. Its coordinate system must be the same of \code{river}. Consider using \code{terra::project} to change
coordinate system.} 
\item{river}{River object. It must have been aggregated (via \code{aggregate_river}).}
\item{categorical}{Logical. Is the covariate categorical (e.g. land cover classes)? If \code{x} consists of several
layers, it is possible to specify \code{categorical} as a vector of logicals (one per each layer); alternatively, 
the single value of \code{categorical} is applied to all layers of \code{x}.}
}

\value{
A river object. The following elements are added:
\item{SC$locCov}{Data frame of covariates evaluated as local values for each subcatchment (i.e., mean covariate value
within a catchment).}
\item{SC$upsCov}{Data frame of covariates evaluated as upstream-averaged values for each subcatchment (i.e., mean 
covariate value within the area upstream of a given subcatchment, including the subcatchment itself).}  
}

\details{
If \code{categorical = TRUE}, the number of columns of \code{SC$locCov}, \code{SC$upsCov} is equal to the number of 
unique values of \code{x} within the catchment. Column names are composed as \code{"y_z"}, where  \code{y = names(x)} and 
\code{z} are the unique values of \code{x}. Values correspond to the fraction of pixels (FD nodes) within the local/upstream
area that are covered by a given category (e.g., land cover type).

If \code{categorical = FALSE}, \code{SC$locCov} and \code{SC$upsCov} have a single column named \code{names(x)}. Values 
correspond to the mean covariate value within the local/upstream reference area.

If \code{x} has multiple layers, columns in the data frames are added sequentially. The same occurs if \code{covariate_river}
is run repeated times (for instance, to compute covariates for one SpatRaster object at a time).
}

\seealso{
\code{\link{rast}}, \code{\link{project}}
}

\examples{
\dontshow{if (interactive() && traudem::can_register_taudem()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\donttest{
 
fp <- system.file("extdata/wigger.tif", package="rivnet")
river <- extract_river(outlet=c(637478,237413),
                       DEM=fp)
river <- aggregate_river(river)

# land cover raster file (categorical)
r1 <- terra::rast(system.file("extdata/landcover.tif", package="rivnet")) 
# legend: 1-urban; 2-agriculture; 3-forest; 4-improductive

river <- covariate_river(r1, river)

plot(river$SC$locCov[ , 1], river) # fraction of urban area within a subcatchment
plot(river$SC$upsCov[ , 1], river) # fraction of upstream-averaged urban area  
 
# mean air temperature raster file (continuous) 
r2 <- terra::rast(system.file("extdata/temperature.tif", package="rivnet"))

river <- covariate_river(r2, river, categorical = FALSE)

plot(river$SC$locCov[, 5], river) # the layer has been added after the 4 previous ones
 
names(river$SC$locCov) 
 
 }
\dontshow{\}) # examplesIf}	
 
}
