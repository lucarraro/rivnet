\name{locate_site}
\alias{locate_site}
\title{Locate site in a river}
\usage{
locate_site <- function(X,Y,river, eucledian = TRUE, showPlot = FALSE, 
						xlim = NULL, ylim = NULL)
}
\description{
Finds location of a site (with coordinates X,Y) within a river object.
}

\arguments{
\item{X, Y}{Longitudinal and latitudinal coordinates of the site. They must be given in the same
coordinate system as the \code{river} object.}
\item{river} River object where the site is to be located. It must have been aggregated (via \code{aggregate_river}).}
\item{eucledian}{Logical. Should the location of the site in the \code{river} be attributed "as the crow flies"? 
If \code{FALSE}, the downstream path from (\code{X}, \code{Y}) is followed until the river network is met.}
\item{showPlot}{Logical. Should a close-up of the relocated site be shown?}
\item{xlim, ylim}{Ranges of x- and y-axis limits for the plot. Only active if \code{showPlots = TRUE}.}
}

\value{
A list with objects:
\item{RNode}{index at the RN level of the relocated site.}
\item{AGode}{index at the AG level of the relocated site.}
}

\details{
This function identifies the node in the river network (at the RN and AG levels) that is closest to
an arbitrary site of coordinates \code{X}, \code{Y}. 

Desired coordinates X, Y can be found in an interactive way by clicking on the river map and using 
function \code{locator}.
}

\seealso{
\code{\link{locator}}
}

\examples{
 fp <- system.file("extdata/wigger.tif", package = "rivers")
 r <- extract_river(outlet = c(637478, 237413),
	                  DEM = fp)					  
 r <- aggregate_river(r)
 

}

}