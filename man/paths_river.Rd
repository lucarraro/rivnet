\name{paths_river}
\alias{paths_river}
\title{Find paths in a river}
\usage{
paths_river(river, ...)
}
\description{
Find paths in a river
}
\arguments{
\item{river}{A \code{river} object.}
\item{...}{Further arguments to be passed to \code{aggregate_OCN}.}
}
\value{
A \code{river} object. See \code{paths_OCN} for description of its structure.}
\details{
This is an alias to \code{\link{paths_OCN}}.
}
\examples{
\dontshow{if (interactive() && traudem::can_register_taudem()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\donttest{
 fp <- system.file("extdata/wigger.tif", package="rivnet")
 r <- extract_river(outlet=c(637478,237413),
	DEM=fp)
r <- aggregate_river(r)	
r <- paths_river(r)
}
\dontshow{\}) # examplesIf}
}