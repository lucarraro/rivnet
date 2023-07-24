\name{points_colorscale}
\alias{points_colorscale}
\title{Draw points with a colorscale}
\usage{
points_colorscale(X, Y, values,
                  bg.palette = hcl.colors(1000, "Reds 3", rev=T),
                  col.palette = hcl.colors(1000, "Reds 3",rev=T),
                  bg.range = NULL, col.range = NULL,
                  pch = 21, cex = 2, lwd = 1.5, force.range = TRUE,
                  add.col.legend = FALSE)
}
\description{
Draw points with values displayed as colors. Two different sets of values can be shown simultaneously: 
one in the background and one in the contour.
}
\arguments{
\item{X, Y}{Coordinates of the points to be displayed.}
\item{values}{values of the quantity to be shown at the sampling sites. It can be a single vector
        (colors are shown on the inside of the points), or a data frame (the first column
         is related to colors on the inside, the second to colors on the outside of the points).}
\item{bg.palette, col.palette}{Color palettes for values on the inside and outside of the points, respectively.}
\item{bg.range, col.range}{Ranges for the legend for values on the inside and outside of the points, 
							respectively. If not specified, the range of values (min-max) is used.}
\item{pch, cex, lwd}{Same as in plot() (Note: only use values between 21-25 for pch).}
\item{force.range}{Locical. If TRUE, values outside the range are constrained at the boundaries of 
the range; if FALSE, a transparent color is used.}
\item{add.col.legend}{Logical. Add a legend for values on the outside of the points?}							
}
\value{
No output is produced. This function is used for its side effetcs.
}
\details{
A call to \code{points} is performed in the background. Therefore, a plot window must be open 
when this function is called.
}
\examples{
fp <- system.file("extdata/wigger.tif", package="rivnet")
river <- extract_river(outlet=c(637478,237413), DEM=fp)
river <- aggregate_river(river)

# some random location of sampling sites, just to test the function
samplingSites <- c(2,15,30,78,97,117,132,106,138,153,156,159,
                    263,176,215,189,11,70,79,87,45,209,26,213)
# we use drainage area as an example variable to be shown

# 1) the function must be called after "plot(river)"
plot(river)
points_colorscale(river$AG$X[samplingSites], river$AG$Y[samplingSites],
                   river$AG$A[samplingSites])

# 2) change color palette
plot(river)
points_colorscale(river$AG$X[samplingSites], river$AG$Y[samplingSites],
                  river$AG$A[samplingSites],
                  bg.palette = hcl.colors(1000, "Inferno"))

# 3) impose a different range
plot(river)
points_colorscale(river$AG$X[samplingSites], river$AG$Y[samplingSites],
                  river$AG$A[samplingSites],
                  bg.range = c(0, 1e8))

# 4) show values outside the range as transparent
plot(river)
points_colorscale(river$AG$X[samplingSites], river$AG$Y[samplingSites],
                  river$AG$A[samplingSites],
                  bg.range = c(0, 1e8), force.range = FALSE)

# 5) show values both on inside and outside of the points (
# drainage area at the upstream vs. downstream end of the reach)
plot(river)
points_colorscale(river$AG$X[samplingSites], river$AG$Y[samplingSites],
            data.frame(river$AG$A[samplingSites], 1.5*river$AG$A[samplingSites]),
            bg.range = c(0, 1e8), col.range = c(0, 1e8), 
            lwd = 4)# increase contour line so it's more visible
# specify same range for both bg.range and col.range
# otherwise they will be shown on different scale			


# 6) same as before, but show two different quantities: 
# drainage area (inside) vs. elevation (outside)
# use different color palettes and add legend for the second color palette
plot(river)
points_colorscale(river$AG$X[samplingSites], river$AG$Y[samplingSites],
            data.frame(river$AG$A[samplingSites], river$AG$Z[samplingSites]),
            col.palette = terrain.colors(1000),
            lwd = 4, add.col.legend = TRUE)
}
