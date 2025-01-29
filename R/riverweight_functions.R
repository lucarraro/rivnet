rast_riverweight <- function(x, river,
                             categorical = TRUE,
                             weightNum = list(func="exponential",
                                              mode="flow-based",
                                              dist50=500,
                                              stream=FALSE,
                                              FA=FALSE),
                             weightDen = NULL){

  allowed_names <- c("mode",   "stream", "func",   "dist50", "FA", "distLinear", "distCauchy", "distExponential", "distPower")
  if( !all(names(weightNum) %in% allowed_names)) stop('Invalid names in weightNum')
  if( !all(names(weightDen) %in% allowed_names)) stop('Invalid names in weightDen')

  if (is.null(weightNum$func)) weightNum$func <- "exponential"
  if (is.null(weightNum$mode)) weightNum$mode <- "flow-based"
  if (is.null(weightNum$stream)) weightNum$stream <- FALSE
  if (is.null(weightNum$FA)) weightNum$FA <- FALSE
  if (length(c(weightNum$dist50, weightNum$distExponential, weightNum$distLinear,
               weightNum$expPower, weightNum$distCauchy)) == 0) {weightNum$dist50 <- 500}

  for (nam in names(weightNum)){
    if (is.null(weightDen[[nam]])) weightDen[[nam]] <- weightNum[[nam]]
  }

  if (weightNum$stream!=weightDen$stream) stop('weightNum$stream must be equal to weightDen$stream')
  if (weightNum$FA!=weightDen$FA) stop('weightNum$FA must be equal to weightDen$FA')
  if (weightNum$mode!=weightDen$mode) stop('weightNum$mode must be equal to weightDen$mode')

  weightNum <- eval_scale.length(weightNum)
  weightDen <- eval_scale.length(weightDen)

  # create raster on catchment extent
  river_rast <- rast(xmin=river$xllcorner - 0.5*river$cellsize,
                     ymin=river$yllcorner - 0.5*river$cellsize,
                     xmax=river$xllcorner - 0.5*river$cellsize + river$dimX*river$cellsize,
                     ymax=river$yllcorner - 0.5*river$cellsize + river$dimY*river$cellsize,
                     resolution=river$cellsize)

  if (categorical[1]){method <- "near"} else {method <- "bilinear"}

  r3 <- resample(x, river_rast, method=method)
  LC_values <- values(r3)
  nlyrs <- dim(LC_values)[2]
  if (nlyrs > length(categorical)){ categorical <- rep(categorical, nlyrs)}

  nms <- NULL
  val <- NULL

  # create matrix with one column per each layer and category
  for (lyr in nlyrs){
    val_all <- values(r3)[river$FD$toDEM, lyr]
    if (categorical[lyr]){
      cats <- sort(unique(val_all))
      val.tmp <- matrix(0,length(val_all), length(cats))
      for (icat in 1:length(cats)){
        val.tmp[val_all==cats[icat],icat] <- 1
        nms <- cbind(nms, paste0(names(r3)[lyr],"_",as.character(cats[icat])))
      }
      val <- cbind(val, val.tmp)
    } else {val <- cbind(val,val_all)
    nms <- cbind(nms, names(r3)[lyr])
    }
  }

  if ((identical(weightNum[order(names(weightNum))],weightDen[order(names(weightDen))]) |
       weightDen$func=="unweighted") &
      weightNum$func=="exponential" & weightNum$mode=="flow-based"){

    val.wu <- eval_wu_exp_flow(val, river, weightNum, weightDen)
    # what if denominator is unweighted? there could still be the fast option
  } else if (weightNum$mode=="flow-based"){
    val.wu <- eval_wu_generic_flow(val, river, weightNum, weightDen)
  } else if (weightNum$mode=="euclidean"){
    if (weightNum$stream){
      val.wu <- eval_wu_euclidean_stream(val, river, weightNum, weightDen)
    } else {
      val.wu <- eval_wu_euclidean(val, river, weightNum, weightDen)
    }}

  val.wu <- cbind(river$FD$X, river$FD$Y, val.wu)
  nms <- cbind("x","y",nms)
  dd <- data.frame(val.wu)
  names(dd) <- nms

  r <- rast(dd, type="xyz", extent=ext(r3))

  invisible(r)
}

eval_wu_exp_flow <- function(val, river, weightNum, weightDen){
  leng <- river$FD$leng
  if (weightNum$stream){leng[river$RN$toFD] <- 0 }
  wl <- exp(-leng/weightNum$scale.length) # doesn't handle the case where den has different scale length
  if (!is.matrix(val)) val <- matrix(val)
  if (weightDen$func=="unweighted") {unweighted=TRUE} else {unweighted=FALSE}
  val_wu <- eval_wu_exp_cpp(val, river, wl, weightNum$FA, unweighted)
  invisible(val_wu)
}

eval_wu_euclidean <- function(val, river, weightNum, weightDen){
  hw <- which(colSums(river$FD$W)==0)
  if (!is.matrix(val)) val <- matrix(val)
  if (weightDen$func=="unweighted") {unweighted=TRUE} else {unweighted=FALSE}
  if (identical(weightNum[order(names(weightNum))],weightDen[order(names(weightDen))]) | unweighted){
    val_wu <- eval_wu_euclidean_cpp_equalND(val,river, weightNum, hw, unweighted)
  } else {
    val_wu <- eval_wu_euclidean_cpp(val, river,
                                    weightNum, weightDen, hw)
  }

  invisible(val_wu)
}

eval_wu_euclidean_stream <- function(val, river, weightNum, weightDen){
  hw <- which(colSums(river$FD$W)==0)
  no_river <- which(river$FD$toRN==0)
  rvr <- which(river$FD$toRN!=0)
  if (!is.matrix(val)) val <- matrix(val)
  if (weightDen$func=="unweighted") {unweighted=TRUE} else {unweighted=FALSE}
  distRiver <- dist_to_river_cpp(river, no_river, rvr)
  if (identical(weightNum[order(names(weightNum))],weightDen[order(names(weightDen))]) | unweighted){
    val_wu <- eval_wu_euclidean_stream_cpp_equalND(val, river, weightNum,
                                                   hw, distRiver, rvr, unweighted)
  } else {
    val_wu <- eval_wu_euclidean_stream_cpp(val, river, weightNum, weightDen,
                                           hw, distRiver, rvr)
  }
  invisible(val_wu)
}

eval_wu_generic_flow <- function(val, river, weightNum, weightDen){
  hw <- which(colSums(river$FD$W)==0)
  if (!is.matrix(val)) val <- matrix(val)
  if (weightDen$func=="unweighted") {unweighted=TRUE} else {unweighted=FALSE}
  if (identical(weightNum[order(names(weightNum))],weightDen[order(names(weightDen))]) | unweighted){
  val_wu <- eval_wu_generic_flow_cpp_equalND(val, river, weightNum, hw, unweighted)
  } else {
    val_wu <- eval_wu_generic_flow_cpp(val, river, weightNum, weightDen, hw)
  }
  invisible(val_wu)
}





get_riverweight <- function(x, rst, river, args_locate_site=list()){ # x: point; rst: raster


  if ("RNnode" %in% names(x)){
    # do nothing
    pp <- x$RNnode
  } else if (is.matrix(x) || is.data.frame(x)){# df or matrix case
    if (ncol(x)==2){
      pp <- numeric(nrow(x))
      for (i in 1:nrow(x)){
        args_locate_site$X <- x[i,1]
        args_locate_site$Y <- x[i,2]
        args_locate_site$river <- river
        tmp <- do.call(locate_site,args_locate_site)
        #tmp <- locate_site(x[i,1],x[i,2],river)
        pp[i] <- tmp$RNnode
      }
    } else { stop("When x is a matrix or data frame, it must have 2 columns.")}
  } else if (is.numeric(x)){# vector case
    if (length(x)==2){
      args_locate_site$X <- x[1]
      args_locate_site$Y <- x[2]
      args_locate_site$river <- river
      tmp <- do.call(locate_site,args_locate_site)
      #tmp <- locate_site(x[1],x[2],river)
      pp <- tmp$RNnode
    } else {stop("When x is a numeric, it must have a length of 2.")}
  } else {
    pp <- numeric(length(x))
    for (i in 1:length(x)){ # this behaves different for matrices vs. vectors vs, data.frames
      tmp <- x[[i]]
      if (("RNnode" %in% tmp)){
        pp[i] <- tmp$RNnode
        # call internal function that does the job? otherwise pp is rewritten
      } else {stop("Invalid format for x.")}
    }
  }

  X_FD <- river$FD$X[river$RN$toFD[pp]]
  Y_FD <- river$FD$Y[river$RN$toFD[pp]]
  cc <- cellFromXY(rst,matrix(c(X_FD,Y_FD),length(pp),2))

  ll <- rst[cc]
  invisible(ll)
}


eval_scale.length <- function(weightList){
  if (length(c(weightList$dist50, weightList$distExponential, weightList$distLinear,
               weightList$expPower, weightList$distCauchy)) > 1){
    stop('Only one parameter among dist50, distExponential, distLinear,
         expPower, distCauchy can be specified')
  }

  if (!is.null(weightList$dist50)){
    if (weightList$func=="exponential" | weightList$func=="gexponential") weightList$scale.length <- weightList$dist50/log(2)
    if (weightList$func=="cauchy") weightList$scale.length <- weightList$dist50
    if (weightList$func=="power") weightList$scale.length <- log(2)/log(1+weightList$dist50) # this is not exactly a length
    if (weightList$func=="linear") weightList$scale.length <- 2*weightList$dist50
  } else {
    if (weightList$func=="exponential"| weightList$func=="gexponential") weightList$scale.length <- weightList$distExponential
    if (weightList$func=="cauchy") weightList$scale.length <- weightList$distCauchy
    if (weightList$func=="power") weightList$scale.length <- weightList$expPower # this is not exactly a length
    if (weightList$func=="linear") weightList$scale.length <- weightList$distLinear
  }
  if (weightList$func=="unweighted") weightList$scale.length <- Inf
  if (is.null(weightList$scale.length)) stop('Mismatch between distFunction and respective parameter chosen.')
  invisible(weightList)
}



