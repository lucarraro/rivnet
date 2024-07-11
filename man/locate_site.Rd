\name{locate_site}
\alias{locate_site}
\title{Locate site in a river}
\usage{
locate_site(X, Y = NULL, river, euclidean = TRUE, showPlot = FALSE, 
  xlim = NULL, ylim = NULL)
}
\description{
Finds location of a site (with coordinates X, Y) within a \code{river} object.
}

\arguments{
\item{X}{Either a list or a numeric value. If \code{X} is a list, it must contain the longitudinal and latitudinal 
coordinates of the site in columns \code{x} (or \code{X}) and \code{y} (or \code{Y}), respectively. If numeric, it represents 
the longitudinal coordinate of the site. Coordinate values must be given in the same coordinate system as the \code{river} object.}
\item{Y}{Latitudinal coordinate of the site. If \code{X} is a list, it can contain the \code{river} object.}
\item{river}{\code{river} object where the site is to be located. It must have been aggregated (via \code{aggregate_river}).}
\item{euclidean}{Logical. Should the location of the site in the \code{river} be attributed "as the crow flies"? 
If \code{FALSE}, the downstream path from (\code{X}, \code{Y}) is followed until the river network is met.}
\item{showPlot}{Logical. Should a close-up of the relocated site be shown?}
\item{xlim, ylim}{Ranges of x- and y-axis limits for the plot. Only active if \code{showPlot = TRUE}.}
}

\value{
A list with objects:
\item{FDode}{index at the FD level of the node that is closest to (\code{X}, \code{Y}). This is generally not a node
that belongs to the river network.}
\item{distance}{The distance between \code{FDnode} and \code{RNnode}, expressed either as Euclidean distance (if 
\code{euclidean = TRUE}), or as downstream distance alternatively.}
\item{RNode}{index at the RN level of the relocated site.}
\item{AGode}{index at the AG level of the relocated site.}
\item{par}{List of graphical parameters as produced by a call to \code{\link{par}} 
(only if \code{showPlot = TRUE}). This can be used to add features to the existing plot. Deprecated.}
}

\details{
This function identifies the node in the river network (at the RN and AG levels) that is closest to
an arbitrary site of coordinates \code{X}, \code{Y}. Only a single site can be processed per function call.

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
 out1 <- locate_site(X, Y, r, showPlot = TRUE) # as the crow flies
 out2 <- locate_site(X, Y, r, showPlot = TRUE, euclidean = FALSE) # follow downstream path
 
 }
 
  # define X, Y by clicking on the map
 if (interactive()) {
 fp <- system.file("extdata/wigger.tif", package = "rivnet")
 r <- extract_river(outlet = c(637478, 237413),
	                  DEM = fp)					  
 r <- aggregate_river(r)
 plot(r)
 
 point <- locator(1) # click on the map to define point
 locate_site(point$X, point$Y, r)
 
 # alternative: specify X as a list and pass river as second argument
 locate_site(point, r)
 }
\dontshow{\}) # examplesIf}
}
