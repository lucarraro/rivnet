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
\item{bg.palette, col.palette}{Color palettes for values on the inside and outside of the points, 
		respectively. Note: Don't change the number of colors (1000); see last example.}
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
