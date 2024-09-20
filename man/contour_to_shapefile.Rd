\name{contour_to_shapefile}
\alias{contour_to_shapefile}
\title{Export catchment contour as shapefile}
\usage{
contour_to_shapefile(river, filename,  
                    EPSG = NULL, ...)
}
\description{
Export catchment contour as shapefile. 
}
\arguments{
\item{river}{\code{river} object.}
\item{filename}{Character. Output file name. It must contain the ".shp" extension.}
 \item{EPSG}{EPSG code. Default is \code{NULL}, which does not produce a .prj file (i.e., the shapefile 
 does not contain projection information).}	
\item{...}{Additional arguments to be passed to \code{writeVector} (e.g., \code{overwrite = TRUE} 
allows overwriting an existing shapefile).}						
}
\value{
No output is produced. This function is used for its side effetcs.
}
\seealso{
\code{\link{writeVector}}
}
\examples{
\dontshow{if (interactive() && traudem::can_register_taudem()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
library(terra) # to use "vect"
fp <- system.file("extdata/wigger.tif", package="rivnet")
river <- extract_river(outlet=c(637478,237413), DEM=fp)

tmpname <- paste0(tempfile(), ".shp")
contour_to_shapefile(river, tmpname, overwrite = TRUE)

# read output
vv <- vect(tmpname)
vv
plot(vv)

\donttest{
# export contour shapefile for multiple catchments
river <- extract_river(outlet=data.frame(x=c(637478,629532),y=c(237413,233782)),
                   EPSG=21781, #CH1903/LV03 coordinate system
                   ext=c(6.2e5,6.6e5,2e5,2.5e5),
                   z=8)
				   
contour_to_shapefile(river, tmpname, overwrite = TRUE)
vv <- vect(tmpname)
vv
plot(vv)

# add projection 
contour_to_shapefile(river, tmpname, 
					EPSG = 21781, 
					overwrite = TRUE)
vv <- vect(tmpname)
vv					
}

\dontshow{\}) # examplesIf}
}
