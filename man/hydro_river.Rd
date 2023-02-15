\name{hydro_river}
\alias{hydro_river}
\title{Assign hydraulic variables to a river network}
\usage{
hydro_river(x, river, level = "AG", leopold = TRUE,
 expWidth = 0.5, expDepth = 0.4, expQ = 1,
 crossSection = "natural", ks = 30, minSlope = NULL)
}
\description{
Assign hydraulic variables (width, water depth, discharge, water velocity, ...) across a river network from measured values based on 
scaling relationships and/or uniform flow equations.
}

\arguments{
\item{x}{Data frame containing measured hydraulic variables. It must consist of columns \code{data} (numeric values of width in m, depth in m or discharge in m^3 s^(-1)), \code{type} (containing the variable type: width (\code{"w"}), depth (\code{"d"}) or discharge (\code{"Q"})), and \code{node} (ID of AG or RN nodes to where \code{data} have been measured). See details.} 
\item{river}{River object. It must have been aggregated (via \code{aggregate_river}).}
\item{level}{Aggregation level at which the nodes in \code{x$node} are defined and at which hydraulic variables are calculated. Possible values are \code{"RN"}, \code{"AG"}. See \code{\link{OCNet}} for details.}
\item{leopold}{Logical. Should scaling relationships of hydraulic variables with drainage area (in the spirit of Leopold and Maddock, 1953) be preferred over uniform flow (Gauchler-Strickler/Manning) relationships? See details.}
\item{expWidth}{Exponent for the scaling relationship of width to drainage area. See details.}
\item{expDepth}{Exponent for the scaling relationship of depth to drainage area. See details.} 
\item{expQ}{Exponent for the scaling relationship of discharge to drainage area. See details.} 
\item{crossSection}{Shape of the river cross-section (constant across the river network). Possible values are \code{"rectangular"}, \code{"natural"}, or a numeric value indicating the exponent of the width-depth relationship. See details.}
\item{ks}{Roughness coefficient according to Gauchler-Strickler (in m^(1/3)s^(-1)). It is the inverse of Manning's roughness coefficient. It can be a single value (thus assumed constant for the whole river network), or a vector of length equal to the number of nodes at the specified \code{level} (one roughness coefficient value for each corresponding node).} 
\item{minSlope}{Minimum slope value, replacing null or \code{NaN} values of slope in the river object. If \code{NULL}, it is assumed equal to either \code{river$slope0} (if the river is an OCN) or to the minimum positive slope value at the selected \code{level}.}
}
\value{
A river object. The following elements are added to the list indicated by \code{level}:
\item{width, depth, discharge}{Values assigned at all nodes (see Details). Units are m, m, m^3 s^(-1), respectively.}
\item{velocity}{Values in m s^(-1). Calculated by continuity (i.e., ratio between discharge and cross-sectional area).} 
\item{volume}{Values in m^3. Calculated as cross-sectional area times length.}
\item{hydraulicRadius}{Values in m. Calculated as the ratio between cross-sectional area and wetted perimeter.}
\item{shearStress}{Values in N m^(-2). Shear stress exerted at the streambed by waterflow. Calculated as gamma*hydraulicRadius*slope, where gamma = 9806 N m^(-3) is the specific weight of water.}
}

\details{
This function is a more complete version of \code{rivergeometry_OCN}.

\code{x} must consist of at least one width value and one value of either depth or discharge. All values included in \code{x$data} must be referred to the same time point, so that spatial interpolation can be performed. If the goal is assessing spatio-temporal changes in hydraulic variables, then \code{hydro_river} must be run independently for each time point. For each node, one cannot specify multiple values of the same variable type. Function \code{locate_site} can be used to attribute \code{x$node}.

If \code{level = "AG"}, the drainage area values used for the power law relationships are calculated as \code{0.5*(river$AG$A + river$AG$AReach)}. If \code{level = "RN"}, \code{river$RN$A} is used.

Width values in \code{x} are assumed to be measured at the water surface. If a single width value is provided, widths are calculated at all nodes from a power-law relationship on drainage area with exponent \code{expWidth} and such that width at the measured node is equal to the provided value. If multiple values of width are provided, the function fits a width-drainage area power law on the provided values. In this case, \code{expWidth} is not used and the output (fitted) width values at the measured nodes are generally different than the observed ones.  

Depending on the type of depth and discharge data in \code{x}, the function behaves in eight different ways:
\enumerate{
\item If one depth value and zero discharge values are provided, the Gauchler-Strickler uniform flow relationship (hereafter GS) is applied to find discharge at the node where depth was measured. Discharge values are then attributed to all nodes based on a power-law relationship vs. drainage area (hereafter PL) with exponent \code{expQ}. Finally, depth values at all nodes are derived from GS.
\item If one discharge value and zero depth values are provided, discharge values are attributed to all nodes based on a PL with exponent \code{expQ}, and such that the value at the measurement node be equal to the observed one. Depth values at all nodes are then derived from GS. 
\item If one discharge and one depth value are provided (not necessarily referred to the same node), discharge values are first attributed as in case 2. If \code{leopold = TRUE}, depth values are derived from a PL with exponent \code{expDepth}; conversely, GS is applied. 
\item If multiple values of discharge and zero values of depth are provided, discharge values are attributed from a power-law fit on measured values vs. drainage area (heareafter PLF). Depth values are then obtained from GS.
\item If multiple values of depth and zero values of discharge are provided, depth values are obtained by PLF. Discharge values are then calculated from GS.
\item If multiple values of discharge and one value of depth are provided, discharge values are first computed as in case 4. Depth values are then obtained by either PF with exponent \code{expDepth} (if \code{leopold = TRUE}), or alternatively via GS.
\item If multiple values of depth and one value of discharge are provided, depth values are first computed as in case 5. Discharge values are then obtained as in case 2 (if \code{leopold = TRUE}), or alternatively computed from GS.  
\item If multiple values of both discharge and depth are provided, discharge values are computed as in case 4, and depth values are computed as in case 5.    
}

Cross-sections are assumed as vertically symmetric. If \code{crossSection = "natural"}, the relationship between width and depth at a cross-section is expressed by width ~ depth^0.65, as suggested by Leopold and Maddock (1953) (where width ~ discharge^0.26 and depth ~ discharge^0.4). Assuming \code{crossSection = 0} is equivalent to \code{"rectangular"} (width does not depend on depth), while \code{crossSection = 1} corresponds to an isosceles triangular cross-section.
}

\seealso{
\code{\link{aggregate_river}}, \code{\link{locate_site}}, \code{\link{rivergeometry_OCN}}. 
}

\examples{
\donttest{
 
fp <- system.file("extdata/wigger.tif", package="rivnet")
river <- extract_river(outlet=c(637478,237413),
                       DEM=fp)
river <- aggregate_river(river)

data <- c(12.8,  6.7,  3.3,  1.1,  9.5,  0.8)
type <- c("w",   "w",  "w",  "w",  "Q",  "d")
node <- c( 46,   109,  181,  145,   46,   46) # assume these have been found via locate_site

x <- data.frame(data=data, type=type, node=node)

river1 <- hydro_river(x, river) # case 3
river2 <- hydro_river(x, river, leopold = FALSE) # case 3 (depth calculated via GS)

plot(0.5*(river1$AG$A + river1$AG$AReach), river1$AG$depth) # Power law with exponent 0.4
plot(0.5*(river2$AG$A + river2$AG$AReach), 
	river2$AG$depth) # Higher depths in reaches with small slope

river3 <- hydro_river(x, river, leopold = FALSE, minSlope = 0.002)
plot(0.5*(river1$AG$A + river1$AG$AReach), river1$AG$depth) # Variability is reduced

river <- hydro_river(x[-5, ], river) # case 1
river <- hydro_river(x[-6, ], river) # case 2

 }
 
 
}
