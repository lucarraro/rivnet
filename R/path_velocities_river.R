path_velocities_river <-
  function(river,level=c("RN","AG"),displayUpdates = FALSE){

    for (str in level){
      if (displayUpdates){
        message(sprintf("%s path velocities... \n",str), appendLF = FALSE)
      }

      if (length(river[[str]]$X)==0){
        stop('Missing fields in river. You should run aggregate_river prior to path_velocities_river.') }

      if (is.null(river[[str]]$downstreamPath)){
        stop('Missing paths in river at the desired aggregation level.
           You should run paths_river prior to path_velocities_river.') }

      if (identical(river[[str]]$downstreamPath, river[[str]]$downstreamPathLength)){
        stop('Missing paths. Set option includePaths = TRUE in paths_river.')}

      if (is.null(river[[str]]$velocity)){
        stop('Missing velocities in river at the desired aggregation level.
           You should run hydro_river or rivergeometry_OCN prior to path_velocities_river.') }

      node <- which(river[[str]]$downNode==river[[str]]$outlet)[1] # a node that is not the outlet
      if (river[[str]]$downstreamPathLength[node, node] == river[[str]]$leng[node]){includeDownstreamNode <- TRUE
      } else {includeDownstreamNode <- FALSE}

      ll <- path_vel_cpp(river, str=str, includeDownstreamNode = includeDownstreamNode)
      river[[str]]$pathVelocities <- spam(ll, river[[str]]$nNodes, river[[str]]$nNodes)
    }
    invisible(river)
  }


# pure R version

# path_velocities_river <- function(river, level = c("RN","AG"),
#                                   displayUpdates = FALSE){
#
#   if (length(river$RN$X)==0){
#     stop('Missing fields in river. You should run aggregate_river prior to path_velocities_river.') }
#
#   if ("RN" %in% level){
#     if (is.null(river$RN$downstreamPath)){
#       stop('Missing paths in river at the desired aggregation level.
#            You should run paths_river prior to path_velocities_river.') }
#
#     if (identical(river$RN$downstreamPath, river$RN$downstreamPathLength)){
#       stop('Missing paths. Set option includePaths = TRUE in paths_river.')}
#
#     if (is.null(river$RN$velocity)){
#       stop('Missing velocities in river at the desired aggregation level.
#            You should run hydro_river or rivergeometry_OCN prior to path_velocities_river.') }
#
#     node <- which(river$RN$downNode==river$RN$outlet)[1] # a node that is not the outlet
#     if (river$RN$downstreamPathLength[node, node] == river$RN$leng[node]){includeDownstreamNode <- TRUE
#     } else {includeDownstreamNode <- FALSE}
#
#     river$RN$pathVelocities <- spam(0,river$RN$nNodes,river$RN$nNodes)
#     set_nodes <- matrix(0,river$RN$nNodes^2,2)
#     set_values <- numeric(river$RN$nNodes^2)
#     k <- 1
#     for (i in 1:river$RN$nNodes){
#       for (j in 1:river$RN$nNodes){
#         path <- river$RN$downstreamPath[[i]][[j]]
#         if (!is.null(path) && !(i == river$RN$outlet && j == river$RN$outlet)){
#           if (includeDownstreamNode){
#           set_values[k] <- river$RN$downstreamPathLength[i,j] / (sum(river$RN$leng[path] / river$RN$velocity[path]))
#           } else {
#             if (length(path)>1){
#             set_values[k] <- river$RN$downstreamPathLength[i,j] /
#               (sum(river$RN$leng[path[1:(length(path)-1)]] / river$RN$velocity[path[1:(length(path)-1)]]))
#             } else {set_values[k] <- 0} # if i=j, we use the patch below to attribute velocities
#           }
#           set_nodes[k,] <- c(i,j)
#           k <- k + 1
#         } else if (i == river$RN$outlet && j == river$RN$outlet) {
#           set_values[k] <- river$RN$velocity[i]
#           set_nodes[k,] <- c(i,j)
#           k <- k + 1
#         }
#       }
#       if (displayUpdates){
#         if ((i %% max(1,round(river$RN$nNodes*0.01)))==0){
#         message(sprintf("RN path velocities... %.1f%%\r",i/(1.001*river$RN$nNodes)*100), appendLF = FALSE)}}
#     }
#     set_values <- set_values[-(k:river$RN$nNodes^2)]
#     set_nodes <- set_nodes[-(k:river$RN$nNodes^2),]
#     river$RN$pathVelocities[set_nodes] <- set_values
#     for (i in 1:river$RN$nNodes){
#       river$RN$pathVelocities[i,i] <- river$RN$velocity[i] # patch to correct when length of path is null
#     }
#     if (displayUpdates){
#       message("RN path velocities... 100.0%\n", appendLF = FALSE)}
#   }
#
#   if ("AG" %in% level){
#     if (is.null(river$AG$downstreamPath)){
#       stop('Missing paths in river at the desired aggregation level.
#            You should run paths_river prior to path_velocities_river.') }
#
#     if (identical(river$AG$downstreamPath, river$AG$downstreamPathLength)){
#       stop('Missing paths. Set option includePaths = TRUE in paths_river.')}
#
#     if (is.null(river$AG$velocity)){
#       stop('Missing velocities in river at the desired aggregation level.
#            You should run hydro_river or rivergeometry_OCN prior to path_velocities_river.') }
#
#     node <- which(river$AG$downNode==river$AG$outlet)[1] # a node that is not the outlet
#     if (river$AG$downstreamPathLength[node, node] == river$AG$leng[node]){includeDownstreamNode <- TRUE
#     } else {includeDownstreamNode <- FALSE}
#
#     river$AG$pathVelocities <- spam(0,river$AG$nNodes,river$AG$nNodes)
#     set_nodes <- matrix(0,river$AG$nNodes^2,2)
#     set_values <- numeric(river$AG$nNodes^2)
#     k <- 1
#     for (i in 1:river$AG$nNodes){
#       for (j in 1:river$AG$nNodes){
#         path <- river$AG$downstreamPath[[i]][[j]]
#         if (!is.null(path) && !(i == river$AG$outlet && j == river$AG$outlet)){
#           if (includeDownstreamNode){
#           set_values[k] <- river$AG$downstreamPathLength[i,j] / (sum(river$AG$leng[path] / river$AG$velocity[path]))
#           } else {
#             if (length(path)>1){
#               set_values[k] <- river$AG$downstreamPathLength[i,j] /
#                 (sum(river$AG$leng[path[1:(length(path)-1)]] / river$AG$velocity[path[1:(length(path)-1)]]))
#             } else {set_values[k] <- 0} # if i=j, we use the patch below to attribute velocities
#           }
#           set_nodes[k,] <- c(i,j)
#           k <- k + 1
#         } else if (i == river$AG$outlet && j == river$AG$outlet) {
#           set_values[k] <- river$AG$velocity[i]
#           set_nodes[k,] <- c(i,j)
#           k <- k + 1
#         }
#       }
#       if (displayUpdates){
#         if ((i %% max(1,round(river$AG$nNodes*0.01)))==0){
#         message(sprintf("AG path velocities... %.1f%%\r",i/(1.001*river$AG$nNodes)*100), appendLF = FALSE)}}
#     }
#     set_values <- set_values[-(k:river$AG$nNodes^2)]
#     set_nodes <- set_nodes[-(k:river$AG$nNodes^2),]
#     river$AG$pathVelocities[set_nodes] <- set_values
#     for (i in 1:river$AG$nNodes){
#       river$AG$pathVelocities[i,i] <- river$AG$velocity[i] # patch to correct when length of path is null
#     }
#     if (displayUpdates){
#       message("AG path velocities... 100.0%\n", appendLF = TRUE)}
#   }
#   invisible(river)
# }
