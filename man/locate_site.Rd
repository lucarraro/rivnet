\name{locate_site}
\alias{locate_site}
\title{Locate site in a river}
\usage{
locate_site(X, Y, river, euclidean = TRUE, showPlot = FALSE, 
  xlim = NULL, ylim = NULL)
}
\description{
Finds location of a site (with coordinates X,Y) within a river object.
}

\arguments{
\item{X, Y}{Longitudinal and latitudinal coordinates of the site. They must be given in the same
coordinate system as the \code{river} object.}
\item{river}{River object where the site is to be located. It must have been aggregated (via \code{aggregate_river}).}
\item{euclidean}{Logical. Should the location of the site in the \code{river} be attributed "as the crow flies"? 
If \code{FALSE}, the downstream path from (\code{X}, \code{Y}) is followed until the river network is met.}
\item{showPlot}{Logical. Should a close-up of the relocated site be shown?}
\item{xlim, ylim}{Ranges of x- and y-axis limits for the plot. Only active if \code{showPlots = TRUE}.}
}

\value{
A list with objects:
\item{FDode}{index at the FD level of the node that is closest to (\code{X}, \code{Y}). This is generally not a node
that belongs to the river network.}
\item{distance}{The distance between \code{FDnode} and \code{RNnode}, expressed either as Euclidean distance (if 
\code{euclidean = TRUE}), or as downstream distance alternatively.}
\item{RNode}{index at the RN level of the relocated site.}
\item{AGode}{index at the AG level of the relocated site.}
}

\details{
This function identifies the node in the river network (at the RN and AG levels) that is closest to
an arbitrary site of coordinates \code{X}, \code{Y}. 

Desired coordinates X, Y can be found in an interactive way by clicking on the river map and using 
function \code{locator}.

Nodes at the RN level thus found can be defined as new breakpoints for reaches (see \code{aggregate_OCN} 
and argument \code{breakpoints}).
}

\seealso{
\code{\link{aggregate_OCN}}, \code{\link{locator}}
}

\examples{
\dontshow{if (interactive() && traudem::can_register_taudem()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\donttest{
 fp <- system.file("extdata/wigger.tif", package = "rivnet")
 r <- extract_river(outlet = c(637478, 237413),
	                  DEM = fp)					  
 r <- aggregate_river(r)
 
 X <- 641329; Y <- 227414 
 locate_site(X, Y, r, showPlot = TRUE) # as the crow flies
 locate_site(X, Y, r, showPlot = TRUE, euclidean = FALSE) # follow downstream path
 }
 
 \dontrun{
 # define X, Y by clicking on the map
 fp <- system.file("extdata/wigger.tif", package = "rivnet")
 r <- extract_river(outlet = c(637478, 237413),
	                  DEM = fp)					  
 r <- aggregate_river(r)
 plot(r)
 
 point <- locator(1) # click on the map to define point
 locate_site(point$X, point$Y, r)
 }
\dontshow{\}) # examplesIf}
}
