\name{plot}

\docType{methods}

\alias{plot}
\alias{plot,river,numeric-method}
\alias{plot,numeric,river-method}
\alias{plot,river,missing-method}

\title{Plot a river }

\description{Plots a river object}

\usage{
\S4method{plot}{river,numeric}(x, y, ...)
\S4method{plot}{numeric,river}(y, x, ...)
\S4method{plot}{river,missing}(x, y, ...)
}

\arguments{
\item{x}{A river object.}
\item{y}{A numeric vector to be displayed. It is equivalent to \code{theme} in \code{draw_thematic_OCN}.
If \code{y} is specified, the river must have been aggregated. See details.}
\item{...}{Arguments passed to the plotting functions \code{draw_simple_OCN}, \code{draw_contour_OCN}, 
\code{draw_thematic_OCN}. See details.} 
}

\details{
This is an interface to the plotting functions \code{draw_simple_OCN}, \code{draw_contour_OCN}, 
\code{draw_thematic_OCN}. If the river object does not have an elevation field (i.e., it has been generated
by \code{create_OCN} or \code{create_general_contour_OCN}, but \code{landscape_OCN} has not 
been run), the plotting function used is \code{draw_simple_OCN}. If the elevation field is present, but the river 
has not been aggregated (via \code{aggregate_OCN} or \code{aggregate_river}), the plotting function used is 
 \code{draw_contour_OCN}. If the river has been aggregated,  \code{draw_thematic_OCN} is used.
}

\seealso{
\code{\link{draw_simple_OCN}}, \code{\link{draw_contour_OCN}}, \code{\link{draw_thematic_OCN}}
}

\examples{
fp <- system.file("extdata/wigger.tif", package="rivers")
r <- extract_river(outlet=c(637478,237413),
	DEM=fp)
plot(r)	# equivalent to draw_contour_OCN

\donttest{
r <- aggregate_river(r)
plot(r) # equivalent to draw_thematic_OCN

# equivalent to draw_thematic_OCN with 'theme' specified
plot(r, r$AG$streamOrder, discreteLevels=TRUE) 
}
	
}