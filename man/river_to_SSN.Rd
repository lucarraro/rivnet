\name{river_to_SSN}
\alias{river_to_SSN}
\title{Transform river into SSN object}
\usage{
river_to_SSN(river, ...)
}
\description{
Function that transforms a river object into a SpatialStreamNetwork object.
}
\arguments{
\item{river}{A river object.}
\item{...}{Further arguments to be passed to \code{OCN_to_SSN}.}
}
\details{
This is an alias to \code{\link{OCN_to_SSN}}.
}
\examples{
\dontshow{if (interactive() && traudem::can_register_taudem()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\donttest{
 fp <- system.file("extdata/wigger.tif", package="rivnet")
 r <- extract_river(outlet=c(637478,237413),
	DEM=fp)
r <- aggregate_river(r)	
s <- river_to_SSN(r, level = "AG", obsSites = sample(r$AG$nNodes, 10),
	path = paste(tempdir(),"/river.ssn", sep = ""), importToR = TRUE)
plot(s)
}
\dontshow{\}) # examplesIf}	
}