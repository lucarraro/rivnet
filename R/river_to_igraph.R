river_to_igraph <- function(river, ...){
  g <- OCNet::OCN_to_igraph(river, ...)
  invisible(g)
}
