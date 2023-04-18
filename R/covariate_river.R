covariate_river <- function(x, river,
                            categorical = TRUE,
                            overwrite = FALSE){

  if (length(river$RN$X)==0){
    stop('Missing fields in river. You should run aggregate_river prior to covariate_river.')
  }

  river_rast <- rast(xmin=river$xllcorner - 0.5*river$cellsize,
                     ymin=river$yllcorner - 0.5*river$cellsize,
                     xmax=river$xllcorner - 0.5*river$cellsize + river$dimX*river$cellsize,
                     ymax=river$yllcorner - 0.5*river$cellsize + river$dimY*river$cellsize,
                     resolution=river$cellsize)

  if (categorical[1]){method <- "near"} else {method <- "bilinear"}

  r3 <- resample(x,river_rast,method=method)
  LC_values <- values(r3)

  nlyrs <- dim(LC_values)[2]

  if (nlyrs > length(categorical)){ categorical <- rep(categorical, nlyrs)}

  for (val in 1:nlyrs){
    if (categorical[val]){
      classes <- sort(unique(values(r3)))
      # local covariates
      covariates <- matrix(0,river$SC$nNodes,length(classes))
      for (ind in 1:length(classes)){
        for (i in 1:river$SC$nNodes){
          num <- sum(LC_values[river$FD$toDEM[river$SC$toFD[[i]]]]==classes[ind], na.rm=TRUE)
          den <- sum(!is.na(LC_values[river$FD$toDEM[river$SC$toFD[[i]]]]))
          covariates[i,ind] <- num/den
        }
      }
      locCov <- data.frame(covariates)
      names(locCov) <- paste0(names(r3)[val],"_",as.character(classes))
    } else {
      covariate <- numeric(river$SC$nNodes)
      for (i in 1:river$SC$nNodes){
        covariate[i] <- mean(LC_values[river$FD$toDEM[river$SC$toFD[[i]]]], na.rm=TRUE)
      }
      locCov <- data.frame(covariate)
      names(locCov) <- names(r3)[val]
    }

    # ups covariates
    covariatesUps <- matrix(0,river$SC$nNodes,dim(locCov)[2])
    for (ind in 1:dim(locCov)[2]){
      for (i in 1:river$SC$nNodes){
        tmp <- river$AG$upstream[[i]]
        covariatesUps[i,ind] <- as.numeric(locCov[tmp,ind] %*% river$SC$A[tmp])/sum(river$SC$A[tmp])
      }
    }
    upsCov <- data.frame(covariatesUps)
    names(upsCov) <- names(locCov)

    if (is.null(river$SC$locCov) | overwrite){
      river$SC[["locCov"]] <- locCov
    } else {
      tmp <- river$SC[["locCov"]]
      df <- data.frame(tmp, locCov, check.names=FALSE)
      names(df)<- c(names(tmp),names(locCov))
      river$SC$locCov <- df
    }

    if (is.null(river$SC$upsCov) | overwrite){
      river$SC[["upsCov"]] <- upsCov
    } else {
      tmp <- river$SC[["upsCov"]]
      df <- data.frame(tmp, upsCov, check.names=FALSE)
      names(df)<- c(names(tmp),names(upsCov))
      river$SC$upsCov <- df
    }
  }
  invisible(river)
}
