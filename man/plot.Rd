\name{plot}

\docType{methods}

\alias{plot}
\alias{plot,river,numeric-method}
\alias{plot,numeric,river-method}
\alias{plot,river,missing-method}

\title{Plot a river }

\description{Plots a river object}

\usage{
\S4method{plot}{river,numeric}(x, y, type, ...)
\S4method{plot}{numeric,river}(x, y, type, ...)
\S4method{plot}{river,missing}(x, type, ...)
}

\arguments{
\item{x}{A river object (or a numeric vector if \code{y} is a river).}
\item{y}{A numeric vector to be displayed (or a river if \code{x} is a numeric vector). 
It is equivalent to \code{theme} in \code{draw_subcatchments_OCN} and \code{draw_thematic_OCN}. If \code{y} is specified, the river must 
have been aggregated. See details.}
\item{type}{Optional argument. If \code{type = "SC"} or \code{type = "subcatchments"}, \code{draw_subcatchments_OCN} is used (provided that 
the river object is aggregated); if  \code{type = "elev2D"}, \code{draw_elev2D_OCN} is used; alternatively,
\code{draw_thematic_OCN} is used.}
\item{...}{Arguments passed to the plotting functions \code{draw_simple_OCN}, \code{draw_contour_OCN}, 
\code{draw_thematic_OCN}. See details.} 
}

\details{
This is an interface to the plotting functions \code{draw_simple_OCN}, \code{\link{draw_elev2D_OCN}}, \code{draw_contour_OCN}, \code{draw_subcatchments_OCN},
\code{draw_thematic_OCN}. If the river object does not have an elevation field (i.e., it has been generated
by \code{create_OCN} or \code{create_general_contour_OCN}, but \code{landscape_OCN} has not 
been run), the plotting function used is \code{draw_simple_OCN}. If the elevation field is present, but the river 
has not been aggregated (via \code{aggregate_OCN} or \code{aggregate_river}), the default plotting function used is 
 \code{draw_contour_OCN}. If the river has been aggregated, \code{draw_subcatchments_OCN} or \code{draw_thematic_OCN} are used depending on \code{type}.
 Elevation maps can be produced with \code{type = "elev2D"}, regardless of whether the river has been aggregated.
}

\seealso{
\code{\link{draw_simple_OCN}}, \code{\link{draw_elev2D_OCN}}, \code{\link{draw_contour_OCN}}, \code{\link{draw_subcatchments_OCN}}, \code{\link{draw_thematic_OCN}}
}

\examples{
\dontshow{if (interactive() && traudem::can_register_taudem()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
fp <- system.file("extdata/wigger.tif", package="rivnet")
r <- extract_river(outlet=c(637478,237413),
	DEM=fp)
plot(r)	# equivalent to draw_contour_OCN

\donttest{
r <- aggregate_river(r)
plot(r) # equivalent to draw_thematic_OCN
plot(r,  type = "SC") # equivalent to draw_subcatchments_OCN

# equivalent to draw_thematic_OCN with 'theme' specified
plot(r, r$AG$streamOrder, discreteLevels = TRUE)
plot(r$AG$streamOrder, r, discreteLevels = TRUE)  # swapping arguments is allowed

# equivalent to draw_subcatchments_OCN with 'theme' specified
plot(r, r$SC$Y, type = "SC", addLegend = FALSE)
plot(r$SC$Y, r,  type = "subcatchments", addLegend = FALSE)  # swapping arguments is allowed

# plot elevation map
plot(r, type = "elev2D", drawRiver = TRUE)
}
\dontshow{\}) # examplesIf}	
}