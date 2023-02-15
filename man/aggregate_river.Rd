\name{aggregate_river}
\alias{aggregate_river}
\title{Aggregate a river}
\usage{
aggregate_river(river, ...)
}
\description{
Aggregates a river
}
\arguments{
\item{river}{A river object.}
\item{...}{Further arguments to be passed to \code{aggregate_OCN}.}
}
\details{
This is an alias to \code{\link{aggregate_OCN}}.
}
\examples{
\dontshow{if (interactive() && traudem::can_register_taudem()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\donttest{
 fp <- system.file("extdata/wigger.tif", package="rivnet")
 r <- extract_river(outlet=c(637478,237413),
	DEM=fp)
r <- aggregate_river(r)	
}
\dontshow{\}) # examplesIf}	
}