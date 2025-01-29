\docType{data}
\name{wigger}
\alias{wigger}
\title{River Wigger}
\format{A \code{river} object. See \code{\link{extract_river}} documentation for details.}
\usage{
data(wigger)
}
\description{
It is built via 

\code{wigger <- extract_river(outlet=c(637478,237413),
                    EPSG=21781, 
                    ext=c(6.2e5,6.6e5,2e5,2.5e5),
                    z=9)}

\code{wigger <- aggregate_river(wigger, maxReachLength = 2500)}

\code{hydrodata <- data.frame(data=c(8, 15), type=c("w","Q"), node=wigger$AG$outlet*c(1,1))}

\code{wigger <- hydro_river(hydrodata, wigger)}

\code{r1 <- rast(system.file("extdata/landcover.tif", package="rivnet"))}

\code{wigger <- covariate_river(r1, wigger)}
}
\keyword{datasets}
