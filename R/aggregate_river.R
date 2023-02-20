aggregate_river <- function(river, ...){
  river <- OCNet::aggregate_OCN(river, ...)
  invisible(river)
}
