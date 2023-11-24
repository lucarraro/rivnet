points_colorscale <- function(X, Y, values,
                              bg.palette = hcl.colors(1000, "Reds 3", rev=T),
                              col.palette = hcl.colors(1000, "Reds 3",rev=T),
                              bg.range = NULL, col.range = NULL,
                              pch = 21, cex = 2, lwd = 1.5, force.range = TRUE,
                              add.col.legend = FALSE, ...){

  if (is.data.frame(values)){
    values.bg <- values[,1]
    values.col <- values[,2]
    if (is.null(col.range)) col.range <- range(values.col, na.rm=T)
  } else {
    values.bg <- values
  }
  if (is.null(bg.range)) bg.range <- range(values.bg, na.rm=T)

  bg.length <- length(bg.palette)
  col.length <- length(col.palette)

  for (ind in 1:length(X)){
    ind.bg <- ceiling(bg.length*(values.bg[ind] - bg.range[1])/(bg.range[2] - bg.range[1]))
    if (is.na(ind.bg) | ((ind.bg < 1 | ind.bg > bg.length) & !force.range)){
      bg <- "#00000000"
    } else {
      ind.bg[ind.bg<1] <- 1; ind.bg[ind.bg>bg.length] <- bg.length
      bg <- bg.palette[ind.bg]
    }

    if (is.data.frame(values)){
      ind.col <- ceiling(col.length*(values.col[ind] - col.range[1])/(col.range[2] - col.range[1]))
      if (is.na(ind.col) | ((ind.col < 1 | ind.col > col.length) & !force.range)){
        col <- "#00000000"
      } else {
        ind.col[ind.col<1] <- 1; ind.col[ind.col>col.length] <- col.length
        col <- col.palette[ind.col]
      }
      points(X[ind], Y[ind], bg=bg, col=col, lwd=lwd, cex=cex, pch=pch)
    } else {
      points(X[ind], Y[ind], bg=bg, lwd=lwd, cex=cex, pch=pch)}
  }

  # add legend
  imagePlot(col=bg.palette, add=T,legend.only=T,
                    zlim=bg.range, ...)

  if (add.col.legend & is.data.frame(values)){
    imagePlot(col=col.palette, add=T,legend.only=T,
                      zlim=col.range, horizontal=T, ...)}


  invisible()

}

