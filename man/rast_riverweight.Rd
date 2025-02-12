\name{rast_riverweight}
\alias{rast_riverweight}

\title{Compute upstream-averaged raster}

\description{
Compute a raster file where each cell value is a weighted average of upstream values.
}

\usage{
rast_riverweight(x, river, 
                categorical = TRUE,
                weightNum = list(func = "exponential",
                                 mode = "flow-based",
                                 dist50 = 500,
                                 stream = FALSE,
                                 FA = FALSE),
                weightDen = NULL)
}

\arguments{
  \item{x}{\code{SpatRaster} object (obtained via \code{terra::rast}) containing categorical or continuous variables from which upstream averages are computed. Its coordinate system must be the same of \code{river}. Consider using \code{terra::project} to change
coordinate system.} 
  \item{river}{A \code{river} object generated via \code{aggregate_river}.}
  \item{categorical}{Logical. Is the covariate categorical (e.g. land cover classes)? If \code{x} consists of several
layers, it is possible to specify \code{categorical} as a vector of logicals (one per each layer); alternatively, 
the single value of \code{categorical} is applied to all layers of \code{x}.}
  \item{weightNum}{List defining attributes of the numerator of the weighted average. See details. }
  \item{weightDen}{List defining attributes of the denominator of the weighted average. If \code{NULL}, it is taken equal to \code{weightNum}. See details. }
}

\details{
Lists \code{weightNum} and \code{weightDen} can include arguments \code{func}, \code{mode}, \code{stream}, \code{FA} and one between \code{dist50}, \code{distExponential}, \code{distCauchy}, \code{distLinear}, \code{expPower}. If not all of these arguments are provided, default values for \code{weightNum} are used (see examples).
\describe{
\item{\code{func}}{expresses the type of distance decay function used. It must be equal to one among \code{"exponential"}, \code{"cauchy"}, \code{"linear"}, \code{"power"}. Only for \code{weightDen}, the value \code{"unweighted"} is also allowed. Distance decay functions are defined as follows:
\describe{
\item{\code{"exponential"}}{\eqn{w(d)=\exp(1-d/d_E)}}
\item{\code{"cauchy"}}{\eqn{w(d)=d_C^2/(d^2 + d_C^2)}}
\item{\code{"linear"}}{\eqn{w(d)=\max(1-d/d_L, 0)}}
\item{\code{"power"}}{\eqn{w(d)=1/(1+d)^{e_P}}}
\item{\code{"unweighted"}}{\eqn{w(d)=1}}
}
where \eqn{w} is the weight of a given source cell, \eqn{d} the distance (see \code{mode}) from the source to the target cell, \eqn{d_E}, \eqn{d_C}, \eqn{d_L} and \eqn{e_P} are parameters. 
}
\item{\code{mode}}{expresses the way upstream distances are computed. It must be equal to one between \code{"flow-based"} (distances computed along steepest descent paths) and \code{"euclidean"} (i.e., distances as the crow flies).}
\item{\code{dist50}, \code{distExponential}, \code{distCauchy}, \code{distLinear}, \code{expPower}}{Parameters for the distance decay function expressed in \code{func}. Parameter \code{dist50} is the distance at which \eqn{w = 0.5}, and it can be expressed for any choice of \code{func}. The other parameters are specific to a given type of \code{func}, and are equal to the respective parameters in the formulas above (i.e., \code{distExponential} = \eqn{d_E}, \code{distCauchy} = \eqn{d_C}, \code{distLinear} = \eqn{d_L}, \code{expPower} = \eqn{e_P}). All parameters but \code{expPower} are distances expressed in the same unit as \code{x} and \code{river}. \code{expPower} is a positive, dimensionless value; note that the value of \code{expPower} depends on the unit of \code{x} and \code{river} (e.g., if distances in \code{river} are expressed in km, the same \code{expPower} will yield a different distance decay function than if distances in \code{river} are in m).}
\item{\code{stream}}{Logical. If \code{TRUE}, distances along the river network are not accounted for, that is, only distances (either along the steepest descent path or as the crow flies, depending on \code{mode}) from the source cell to the river network are considered. If \code{mode = "euclidean"}, this corresponds to the shortest planar distance between the source cell and any river network cell. This implies \eqn{d = 0} for all source cells lying in the river network.}
\item{\code{FA}}{Logical. Should flow-contributing areas (expressed as numbers of cells upstream of a source cell--including the source cell itself) be included as a multiplicative factor to \eqn{w}?}
}

To ensure computational efficiency for large and highly resolved rasters (note: it is the cell size of the \code{river} object that matters, not the resolution of \code{x}), it is recommended to use \code{func = "exponential"} (and/or \code{func = "exponential"} for \code{weightDen}) and \code{mode = "flow-based"}. Values of \code{stream} and \code{FA} do not affect computational speed.  
}

\value{
A \code{SpatRaster} object containing as many layers as the number of layers in \code{x}, each possibly multiplied by the number of unique categories featured in the layer (if the layer is categorical). If layer \code{y}  \code{x} is categorical, the corresponding layers in the output \code{SpatRaster} object are named \code{y_z}, where \code{z} is the value of a unique category in the original layer \code{y}. If layer \code{y} in \code{x} is continuous, the corresponding layer in the output \code{SpatRaster} object is also named \code{y}.

The output \code{SpatRaster} object has the same extent and resolution (i.e., cell size) of the input \code{river} object.  
}

\seealso{
\code{\link{aggregate_river}}, \code{\link[terra:rast]{terra::rast}}, \code{\link[terra:project]{terra::project}}, \code{\link{get_riverweight}}
}

\examples{
library(terra)
data(wigger)

r1 <- rast(system.file("extdata/landcover.tif", package = "rivnet")) 
# legend: 1-urban; 2-agriculture; 3-forest; 4-improductive

r.exp <- rast_riverweight(r1, wigger)
plot(r.exp)

\donttest{
# unweighted denominator
r.unweighted <- rast_riverweight(r1, wigger, 
					weightDen = list(func = "unweighted"))

# alternative distance decay functions (with same dist50)
# these take more time than the default func = "exponential"
r.cau <- rast_riverweight(r1, wigger, 
						  weightNum = list(func = "cauchy"))
r.lin <- rast_riverweight(r1, wigger, 
						  weightNum = list(func = "linear"))
r.pow <- rast_riverweight(r1, wigger, 
						  weightNum = list(func = "power"))

# ignore distances on the river network
r.exp_S <- rast_riverweight(r1, wigger, 
							weightNum = list(stream = TRUE))

# include flow accumulation in the weight
r.exp_FA <- rast_riverweight(r1, wigger, 
							 weightNum = list(FA = TRUE))

# use Euclidean distances (takes more time)
# Euclidean distance from source to target
r.dO <- rast_riverweight(r1, wigger, 
						 weightNum = list(mode = "euclidean"))
# Euclidean distance from source to river network
r.dOS <- rast_riverweight(r1, wigger, 
						  weightNum = list(mode = "euclidean", 
										   stream = TRUE))

# specify exponential decay parameter in different ways
r.exp1 <- rast_riverweight(r1, wigger, 
						   weightNum = list(dist50 = 1000*log(2)))
r.exp2 <- rast_riverweight(r1, wigger, 
                           weightNum = list(distExponential = 1000))
identical(r.exp1, r.exp2)


}
}