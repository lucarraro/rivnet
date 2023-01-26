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
\item{river}{A river object.}
\item{...}{Further arguments to be passed to \code{aggregate_OCN}.}
}
\details{
This is an alias to \code{\link{paths_OCN}}.
}
\examples{
\donttest{
 fp <- system.file("extdata/wigger.tif", package="rivnet")
 r <- extract_river(outlet=c(637478,237413),
	DEM=fp)
r <- aggregate_river(r)	
r <- paths_river(r)
}
}