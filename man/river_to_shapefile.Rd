\name{river_to_shapefile}
\alias{river_to_shapefile}
\title{Export river network as shapefile}
\usage{
river_to_shapefile(river, filename, atts = NULL, 
                    EPSG = NULL, ...)
}
\description{
Export river network as shapefile. Reach attributes can be added.
}
\arguments{
\item{river}{\code{river} object. It must have been aggregated (via \code{aggregate_river}).}
\item{filename}{Character. Output file name. It must contain the ".shp" extension.}
\item{atts}{Attributes at AG level that can be exported. This should be a character vector, with entries 
 corresponding to the names of fields within the \code{AG} sub-list of \code{river}. See example.}
 \item{EPSG}{EPSG code. Default is \code{NULL}, which does not produce a .prj file (i.e., the shapefile 
 does not contain projection information).}	
\item{...}{Additional arguments to be passed to \code{writeVector} (e.g., \code{overwrite = TRUE} 
allows overwriting an existing shapefile).}						
}
\value{
No output is produced. This function is used for its side effetcs.
}
\seealso{
\code{\link[terra:writeVector]{terra::writeVector}}
}
\examples{
\dontshow{if (interactive() && traudem::can_register_taudem()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
library(terra) # to use "vect"
fp <- system.file("extdata/wigger.tif", package="rivnet")
river <- extract_river(outlet=c(637478,237413), DEM=fp)
river <- aggregate_river(river)

tmpname <- paste0(tempfile(), ".shp")
river_to_shapefile(river, tmpname, overwrite = TRUE)

# read output
vv <- vect(tmpname)
vv
plot(vv)

\donttest{
# add attributes to shapefile (drainage area, stream order)
river_to_shapefile(river, tmpname,
					atts = c("A", "streamOrder"), # same names as in river$AG
					overwrite = TRUE) 
vv <- vect(tmpname)
vv

# add projection 
river_to_shapefile(river, tmpname, 
					EPSG = 21781, 
					overwrite = TRUE)
vv <- vect(tmpname)
vv					
}

\dontshow{\}) # examplesIf}
}
