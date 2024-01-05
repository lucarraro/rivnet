\name{extract_river}
\alias{extract_river}
\title{Extract a river}
\usage{
extract_river(outlet, EPSG=NULL, ext=NULL, z=NULL, DEM=NULL,
  as.river=TRUE, as.rast=FALSE, filename=NULL, showPlot=FALSE,
  threshold_parameter=1000, n_processes=1, displayUpdates=0, src="aws",
  args_get_elev_raster=list())
}
\description{
Function that extracts a river network from elevation data via TauDEM's D8 flow direction
algorithm. It can return a river object and/or output from TauDEM functions as a raster file.
Elevation data can be either downloaded from the web or provided externally.
}

\arguments{
\item{outlet}{A vector, matrix or data frame expressing the coordinates of the river outlet(s)
(in the coordinate system identified by \code{EPSG}, or the same as in \code{DEM}, if the latter
is provided). If a vector, the odd components identify the longitudinal (x) coordinates, and the
even components the latitudinal (y) coordinates. If a matrix, it should have 2 columns (for x and
y coordinates respectively) and number of rows equal to the number of outlets. If a data frame,
it should have components \code{outlet$x}, \code{outlet$y} identifying the respective coordinates.}
\item{EPSG}{EPSG code identifying the coordinate system used. See \url{https://epsg.org/}. It is
required if \code{DEM} is not specified, and not used otherwise. It is recommended to use projected
coordinate systems, so that lengths and areas in the river object are in metric units.}
\item{ext}{Vector expressing the extent of the region where elevation data are downloaded.
It must be in the form \code{c(xmin, xmax, ymin, ymax)}. Coordinates must be expressed in the
coordinate system identified by \code{EPSG}. It is required if \code{DEM} is not specified, and
not used otherwise.}
\item{z}{Zoom level at which elevation data should be downloaded. See \code{get_elev_raster} for
details. It is required if \code{DEM} is not specified, and not used otherwise.}
\item{DEM}{Filename of the Digital Elevation Model raster file to be used.}
\item{as.river}{Logical. Should a river object be created?}
\item{as.rast}{Logical. Should a raster file containing results from \code{traudem} functions output
be created?}
\item{filename}{Filename of the raster file produced if \code{as.rast=TRUE}. Only required if
\code{as.rast=TRUE}. It can be a single filename, or four filenames, in which case four different
raster files are produced. See Details.}
\item{showPlot}{Logical. Should a plot of the calculated contributing area and extracted catchment
contour be produced?}
\item{threshold_parameter}{Value passed to \code{taudem_threshold}. See example and `traudem`
documentation for details.}
\item{n_processes}{Value passed to the \code{traudem} functions.}
\item{displayUpdates}{Numeric. Possible values are \code{0}, \code{1}, \code{2}. If \code{0}, console output 
is suppressed (barring some messages from \code{elevatr::get_elev_raster}, if \code{DEM} is not provided). 
If \code{1}, succint console output is produced. If \code{2}, extensive  console output from 
\code{elevatr::get_elev_raster} and \code{traudem} functions is displayed. }
\item{src}{Value passed to \code{get_elev_raster}. Deprecated.}
\item{args_get_elev_raster}{List of additional parameters to be passed to \code{get_elev_raster}. 
Parameters included in this list override default options (for example, if two different zoom 
values are specified in \code{args_get_elev_raster$z} and \code{z}, then the latter is not considered. )}
}

\value{
A \code{river} object. See \code{create_OCN}, \code{landscape_OCN} for description of its structure.}

\details{
This is a wrapper to \code{elevatr} and \code{traudem} functions, allowing a seamless extraction of
river networks from elevation data. The output \code{river} object is compatible with \code{OCNet} functions
(it is equivalent to an OCN produced by \code{landscape_OCN}).

The workflow of TauDEM commands used is as follows:
 \code{PitRemove} -> \code{D8FlowDir} -> \code{D8ContributingArea} -> \code{StreamDefByThreshold} ->
 \code{MoveOutletsToStreams} -> \code{D8ContributingArea}.
 See \url{https://hydrology.usu.edu/taudem/taudem5/index.html} for details on TauDEM.

 When \code{as.rast = TRUE}, a raster file is returned. It consists of four layers:
 \describe{
 \item{fel}{pit-filled elevation data}
 \item{p}{D8 flow directions}
 \item{ad8}{contributing area for the whole region}
 \item{ssa}{contributing area with respect to the outlet(s) used}
 }
The raster file is written via \code{terra::writeRaster}.

If nested outlets are specified, the function ignores the upstream outlet.
}

\seealso{
\code{\link{get_elev_raster}}, \code{\link{taudem_threshold}},
\code{\link{create_OCN}}, \code{\link{landscape_OCN}}.
}

\examples{
\dontshow{if (interactive() && traudem::can_register_taudem()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# extract the river Wigger (Switzerland) from DEM raster file
# outlet coordinates are expressed in the CH1903/LV03 coordinate system
# (i.e. same as the DEM file)
 fp <- system.file("extdata/wigger.tif", package="rivnet")
 r <- extract_river(outlet=c(637478,237413),
	DEM=fp)
r


\donttest{
# same as above but download DEM data via elevatr
r <- extract_river(outlet=c(637478,237413),
	EPSG=21781, #CH1903/LV03 coordinate system
	ext=c(6.2e5,6.6e5,2e5,2.5e5),
	z=8)
}

\donttest{
# enhance resolution by increasing zoom
r2 <- extract_river(outlet=c(637478,237413),
	EPSG=21781, #CH1903/LV03 coordinate system
	ext=c(6.2e5,6.6e5,2e5,2.5e5),
	z=9)
plot(r)
plot(r2)
}

\donttest{
# specify two outlets as a data frame
r <- extract_river(outlet=data.frame(x=c(637478,629532),y=c(237413,233782)),
                    EPSG=21781, #CH1903/LV03 coordinate system
                    ext=c(6.2e5,6.6e5,2e5,2.5e5),
                    z=10, showPlot=TRUE)
plot(r)

r <- aggregate_river(r)
plot(r, chooseCM = 2)  # display only the second catchment
# (i.e. that identified by the second outlet)
}

\donttest{
# effect of threshold_parameter
r <- extract_river(outlet = c(637478, 237413),
                    EPSG = 21781, #CH1903/LV03 coordinate system
                    ext = c(6.2e5, 6.6e5, 2e5, 2.5e5),
                    z = 8, threshold_parameter = 50,
				    showPlot = TRUE)
plot(r) # if threshold_parameter is too small, the outlet might be located
# in a smaller river reach, and the extracted river network would be too small
# showPlot = TRUE can help identify what is going on

r <- extract_river(outlet = c(637478, 237413),
                    EPSG = 21781, #CH1903/LV03 coordinate system
                    ext = c(6.2e5, 6.6e5, 2e5, 2.5e5),
                    z = 8, threshold_parameter = 1e5,
				    showPlot = TRUE)
plot(r) # if threshold_parameter is too large, the outlet pixel might not be
# located at all (for instance, in this case no cells have contributing area
# above threshold_parameter), hence throwing an error
}
\dontshow{\}) # examplesIf}
}
