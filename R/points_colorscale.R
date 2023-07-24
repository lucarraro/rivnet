points_colorscale <- function(X, Y, values,
                              bg.palette = hcl.colors(1000, "Reds 3", rev=T),
                              col.palette = hcl.colors(1000, "Reds 3",rev=T),
                              bg.range = NULL, col.range = NULL,
                              pch = 21, cex = 2, lwd = 1.5, force.range = TRUE,
                              add.col.legend = FALSE){

  if (is.data.frame(values)){
    values.bg <- values[,1]
    values.col <- values[,2]
    if (is.null(col.range)) col.range <- range(values.col, na.rm=T)
  } else {
    values.bg <- values
  }
  if (is.null(bg.range)) bg.range <- range(values.bg, na.rm=T)

  for (ind in 1:length(X)){
    ind.bg <- ceiling(1000*(values.bg[ind] - bg.range[1])/(bg.range[2] - bg.range[1]))
    if (is.na(ind.bg) | ((ind.bg < 1 | ind.bg > 1000) & !force.range)){
      bg <- "#00000000"
    } else {
      ind.bg[ind.bg<1] <- 1; ind.bg[ind.bg>1000] <- 1000
      bg <- bg.palette[ind.bg]
    }

    if (is.data.frame(values)){
      ind.col <- ceiling(1000*(values.col[ind] - col.range[1])/(col.range[2] - col.range[1]))
      if (is.na(ind.col) | ((ind.col < 1 | ind.col > 1000) & !force.range)){
        col <- "#00000000"
      } else {
        ind.col[ind.col<1] <- 1; ind.col[ind.col>1000] <- 1000
        col <- col.palette[ind.col]
      }
      points(X[ind], Y[ind], bg=bg, col=col, lwd=lwd, cex=cex, pch=pch)
    } else {
      points(X[ind], Y[ind], bg=bg, lwd=lwd, cex=cex, pch=pch)}
  }

  # add legend
  imagePlot(col=bg.palette, add=T,legend.only=T,
                    zlim=bg.range)

  if (add.col.legend & is.data.frame(values)){
    imagePlot(col=col.palette, add=T,legend.only=T,
                      zlim=col.range, horizontal=T)}


  invisible()

}
# X, Y: coordinates of samplign sites
# values: values of the quantity to be shown at the sampling sites. It can be a single vector
#         (colors are shown on the inside of the points), or a data frame (the first column
#         is related to colors on the inside, the second to colors on the outside of the points)
# bg.palette, col.palette: color palettes for values on the inside and outside of the points, respectively
#                          Note: Don't change the number of colors (1000); see last example
# bg.range, col.range: ranges for the legend for values on the inside and outside of the points, respectively.
#                      If not specified, the range of values (min-max) is used
# pch, cex, lwd: as in plot() (pch: points type-only use values between 21-25; cex: points size;
#                lwd: thickness of contour line)
# force.range: if TRUE, values outside the range are constrained at the boundaries of the range;
#               if FALSE, a transparent color is used
# add.col.legend: add a legend for values on the outside of the points?

# # use river from the extract_river example file
# fp <- system.file("extdata/wigger.tif", package="rivnet")
# river <- extract_river(outlet=c(637478,237413), DEM=fp)
# river <- aggregate_river(river)
#
# # some random location of sampling sites, just to test the function
# samplingSites <- c(2,15,30,78,97,117,132,106,138,153,156,159,
#                    263,176,215,189,11,70,79,87,45,209,26,213)
# # we use drainage area as an example variable to be shown
#
# # 1) the function must be called after "plot(river)"
# plot(river)
# points_colorscale(river$AG$X[samplingSites], river$AG$Y[samplingSites],
#                   river$AG$A[samplingSites])
#
# # 2) change color palette
# plot(river)
# points_colorscale(river$AG$X[samplingSites], river$AG$Y[samplingSites],
#                   river$AG$A[samplingSites],
#                   bg.palette=hcl.colors(1000,"Inferno"))
#
#
# # 3) impose a different range
# plot(river)
# points_colorscale(river$AG$X[samplingSites], river$AG$Y[samplingSites],
#                   river$AG$A[samplingSites],
#                   bg.range=c(0, 1e8))
#
# # 4) show values outside the range as transparent
# plot(river)
# points_colorscale(river$AG$X[samplingSites], river$AG$Y[samplingSites],
#                   river$AG$A[samplingSites],
#                   bg.range=c(0, 1e8), force.range=F)
#
# # 5) show values both on inside and outside of the points (drainage area at the upstream vs. downstream end of the reach)
# plot(river)
# points_colorscale(river$AG$X[samplingSites], river$AG$Y[samplingSites],
#                   data.frame(river$AG$A[samplingSites], 1.5*river$AG$A[samplingSites]),
#                   bg.range=c(0, 1e8), col.range=c(0, 1e8), # specify same range for both, otherwise they will be shown on different scales
#                   lwd=4)# increase contour line so it's more visible
#
# # 6) same as before, but show two different quantities: drainage area (inside) vs. elevation (outside)
# # use different color palettes and add legend for the second color palette
# plot(river)
# points_colorscale(river$AG$X[samplingSites], river$AG$Y[samplingSites],
#                   data.frame(river$AG$A[samplingSites], river$AG$Z[samplingSites]),
#                   col.palette = terrain.colors(1000),
#                   lwd=4, add.col.legend=T)
