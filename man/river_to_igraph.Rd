\name{river_to_igraph}
\alias{river_to_igraph}
\title{river_to_igraph}
\usage{
river_to_igraph(river, ...)
}
\description{
Transform a \code{river} in an \code{igraph} object.
}
\arguments{
\item{river}{A \code{river} object.}
\item{...}{Further arguments to be passed to \code{OCN_to_igraph}.}
}
\value{
An \code{igraph} object.}
\details{
This is an alias to \code{\link[OCNet:OCN_to_igraph]{OCNet::OCN_to_igraph}}.
}
\examples{
\dontshow{if (interactive() && traudem::can_register_taudem()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\donttest{
 fp <- system.file("extdata/wigger.tif", package="rivnet")
 r <- extract_river(outlet=c(637478,237413),
	DEM=fp)
r <- aggregate_river(r)	
g <- river_to_igraph(r, level = "AG")
g
}
\dontshow{\}) # examplesIf}	
}