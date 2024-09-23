\name{river_to_AEM}
\alias{river_to_AEM}
\title{river_to_AEM}
\usage{
river_to_AEM(river, ...)
}
\description{
Construct asymmetric eigenvector maps (AEM) from a river
}
\arguments{
\item{river}{A \code{river} object.}
\item{...}{Further arguments to be passed to \code{OCN_to_AEM}.}
}
\value{
A \code{river} object. }
\details{
This is an alias to \code{\link[OCNet:OCN_to_AEM]{OCNet::OCN_to_AEM}}.
}
\examples{
\dontshow{if (interactive() && traudem::can_register_taudem()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\donttest{
 fp <- system.file("extdata/wigger.tif", package="rivnet")
 r <- extract_river(outlet=c(637478, 237413),
	DEM = fp)
r <- aggregate_river(r)
out.aem <- river_to_AEM(r)
}
\dontshow{\}) # examplesIf}
}
