\name{path_velocities_river}
\alias{path_velocities_river}
\title{Calculate velocities along paths in a river}
\usage{
path_velocities_river(river, level = c("RN", "AG"), 
  displayUpdates = FALSE)
}
\description{
Calculate mean water velocities along paths in a river network.
}

\arguments{
\item{river}{River object. It must have been aggregated (via \code{aggregate_river}), and contain paths (via \code{paths_river} with \code{includeDownstreamNode = TRUE}) 
and velocities (via \code{hydro_river} or \code{rivergeometry_OCN}) at the desired aggregation level(s). See \code{level}.}
\item{level}{Aggregation level(s) at which path velocities should be calculated. Possible values are \code{"RN"}, \code{"AG"}, \code{c("RN","AG")}.}
\item{displayUpdates}{Logical. State if updates are printed on the console while \code{path_velocities_river} runs.}
}
\value{
A river object. The following element is added to the list indicated by \code{level}:
\item{pathVelocities}{It is a \code{spam} object. \code{pathVelocities[i, j]} is the mean water velocity along the downstream path from nodes \code{i} to \code{j}, and is null if the two nodes are not connected by such a path.}

}

\details{
Velocities are calculated by dividing the total distance (length of the downstream path joining two nodes) by the total time (sum of times taken to cover all nodes in between the origin and destination nodes; such times are calculated as length/velocity).

Note that paths may or may not include the downstream node; this is controlled by option \code{includeDownstreamNode} in \code{paths_river}. Path velocities are calculated accordingly. 
In both cases, diagonal entries of \code{pathVelocity} are set equal to the respective node velocity. See example.   
}

\seealso{
\code{\link{paths_river}}, \code{\link{hydro_river}}, \code{\link{rivergeometry_OCN}}. 
}

\examples{
\donttest{
 
fp <- system.file("extdata/wigger.tif", package="rivnet")
river <- extract_river(outlet=c(637478,237413),
                       DEM=fp)
river <- aggregate_river(river)
river <- paths_river(river, includePaths = TRUE)
river <- OCNet::rivergeometry_OCN(river) # simplified alternative to hydro_river
                                  # to attribute velocities at all RN and AG nodes
								  
river <- path_velocities_river(river, level = "AG") # downstream nodes are not included in paths
river$AG$pathVelocities[176, 176] 
river$AG$pathVelocities[176, 174]
 # node 174 is immediately downstream of 176; if downstream nodes are not included 
 # in paths, the two velocities are equal

river2 <- paths_river(river, includePaths = TRUE, includeDownstreamNode = TRUE)
river2 <- path_velocities_river(river2, level = "AG") # now downstream nodes are included in paths
river2$AG$pathVelocities[176, 176] 
river2$AG$pathVelocities[176, 174]
 }
 
 
}
