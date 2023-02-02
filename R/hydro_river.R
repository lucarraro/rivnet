hydro_river <- function(x, river, level = "AG",
                        leopold = TRUE,
                        expWidth = 0.5,
                        expDepth = 0.4,
                        expQ = 1,
                        crossSection = "natural",
                        ks = 30,
                        minSlope = NULL){

  if (!("RN" %in% names(river))){
    stop('Missing fields in river. You should run aggregate_river prior to hydro_river.')
  }

  if (!identical(c("data","node","type"), sort(names(x)))){
    stop('Missing or undefined arguments in x.')
  }

  if (typeof(crossSection) == "numeric"){
    alpha <- crossSection
  } else if (crossSection == "rectangular"){
    alpha <- 0
  } else if (crossSection == "natural"){
    alpha <- 0.65
  } else {
    stop('Invalid format for crossSection.')
  }

  if (level =="AG"){
    slope <- river$AG$slope; leng <- river$AG$leng
    A <- 0.5*(river$AG$A+river$AG$AReach); outlet <- river$AG$outlet
  } else if (level=="RN"){
    slope <- river$RN$slope; leng <- river$RN$leng
    A <- river$RN$A; outlet <- river$RN$outlet
  } else {
    stop('Invalid level.')
  }

  # fix outlet variables
  if (is.null(minSlope)){
    if (length(river$slope0)>1){ minSlope <- river$slope0
    } else {
      tmp <- slope[!is.nan(slope)]
      minSlope <- min(slope[slope>0],na.rm=TRUE)}
  }
  leng[leng==0] <- river$cellsize
  slope[is.nan(slope)]  <- minSlope
  slope[slope==0] <- minSlope

  if (length(ks)==1) {ks <- ks+numeric(length(slope))} # one ks value per reach

  data <- x$data
  type <- x$type
  node <- x$node

  width <- data[which(type=="w")]
  node_w <- node[which(type=="w")]
  if (length(unique(node_w)) < length(node_w)){ stop('Only one value of width per node can be provided.')}

  depth <- data[which(type=="d")]
  node_d <- node[which(type=="d")]
  if (length(unique(node_d)) < length(node_d)){ stop('Only one value of depth per node can be provided.')}

  Q <- data[which(type=="Q")]
  node_Q <- node[which(type=="Q")]
  if (length(unique(node_Q)) < length(node_Q)){ stop('Only one value of discharge per node can be provided.')}

  if (length(width)==0){
    stop('At least one width value must be specified.')
  }

  if (length(depth)==0 & length(Q)==0){
    stop('At least one value of either depth or discharge must be specified.')
  }

  # Fit widths
  if (length(width)>1){
    lmod <- summary(lm(log(width) ~ log(A[node_w])))
    width_fitted <- exp(lmod$coefficients[1,1])*A^lmod$coefficients[2,1]
  } else {
    width_fitted <- width*(A/A[node_w])^expWidth
  }

  if (length(width)>1){
    lmod <- summary(lm(log(width) ~ log(A[node_w])))
    width_fitted <- exp(lmod$coefficients[1,1])*A^lmod$coefficients[2,1]
  } else {
    width_fitted <- width*(A/A[node_w])^expWidth
  }

  # WORKFLOW FOR GENERIC CROSS-SECTION
  # width is independent from the rest; then we can have the following cases
  # 1- only one value of depth: apply GS with fixed (or user-attributed) ks -> find velocity and Q in that point
  #   scale Q with power-law (exponent 1 or user-attributed) and determine d from GS relationship -> v from continuity
  # 2- only one value of Q: fit power-law with exponent 1 (or user-attributed) -> use local GS relationships and determine depth
  #    -> determine v by continuity
  # 3- one value of Q and depth at the same (or different) reach(es): fit Leopold power laws with known or user-attributed exponents
  # 4- multiple values of Q: fit lm and derive Q. apply GS to determine depth
  # 5- multiple values of depth: fit lm and derive depth. apply GS to determine Q
  # 6- multiple values of Q and 1 depth: fit lm for Q, GS for depth. velocity by continuity
  #    (why not Leopold for depth? -> depending on a variable type <- c("manning","leopold"))
  # 7- multiple values of d and 1 Q: fit lm for d, GS for Q. velocity by continuity (depending on the variable "type")
  # 8- multiple values of both d, Q: fit lms.

  if (length(depth)==1 & length(Q)==0){  CASE = 1
  Q <- Q_GS_generic(depth, ks[node_d], slope[node_d], width_fitted[node_d], alpha)
  Q_fitted <- Q*(A/A[node_d])^expQ
  depth_fitted <- depth_GS_generic(Q_fitted, ks, slope, width_fitted, alpha)

  } else if (length(depth)==0 & length(Q)==1){  CASE = 2
  Q_fitted <- Q*(A/A[node_Q])^expQ
  depth_fitted <- depth_GS_generic(Q_fitted,ks,slope,width_fitted, alpha)

  } else if (length(depth)==1 & length(Q)==1){  CASE = 3
  Q_fitted <- Q*(A/A[node_Q])^expQ
  if (leopold==TRUE){
    depth_fitted <- depth*(A/A[node_d])^expDepth
  } else {
    depth_fitted <- depth_GS_generic(Q_fitted,ks,slope,width_fitted,alpha)}

  } else if (length(depth)==0){  CASE = 4
  lmod <- summary(lm(log(Q) ~ log(A[node_Q])))
  Q_fitted <- exp(lmod$coefficients[1,1])*A^lmod$coefficients[2,1]
  depth_fitted <- depth_GS_generic(Q_fitted,ks,slope,width_fitted, alpha)

  } else if (length(Q)==0){  CASE = 5
  lmod <- summary(lm(log(depth) ~ log(A[node_d])))
  depth_fitted <- exp(lmod$coefficients[1,1])*A^lmod$coefficients[2,1]
  Q_fitted <- Q_GS_generic(depth_fitted, ks, slope, width_fitted, alpha)

  } else if (length(Q)>1 & length(depth)==1){  CASE = 6
  lmod <- summary(lm(log(Q) ~ log(A[node_Q])))
  Q_fitted <- exp(lmod$coefficients[1,1])*A^lmod$coefficients[2,1]
  if (leopold==TRUE){
    depth_fitted <- depth*(A/A[node_d])^expDepth
  } else {
    depth_fitted <- depth_GS_generic(Q_fitted,ks,slope,width_fitted, alpha)}

  } else if (length(depth)>1 & length(Q)==1){  CASE = 7
  lmod <- summary(lm(log(depth) ~ log(A[node_d])))
  depth_fitted <- exp(lmod$coefficients[1,1])*A^lmod$coefficients[2,1]
  if (leopold==TRUE){
    Q_fitted <- Q*(A/A[node_Q])^expQ
  } else {
    Q_fitted <- Q_GS_generic(depth_fitted, ks, slope, width_fitted, alpha)
  }

  } else if (length(depth)>1 & length(Q)>1){  CASE = 8
  lmod <- summary(lm(log(depth) ~ log(A[node_d])))
  depth_fitted <- exp(lmod$coefficients[1,1])*A^lmod$coefficients[2,1]
  lmod <- summary(lm(log(Q) ~ log(A[node_Q])))
  Q_fitted <- exp(lmod$coefficients[1,1])*A^lmod$coefficients[2,1]

  }
  velocity_fitted <- Q_fitted/depth_fitted/width_fitted*(1+alpha)

  RH <- eval_RH(width_fitted, depth_fitted, alpha)
  tau <- 9806*RH*slope
  volume <- leng*width_fitted/(1+alpha)*depth_fitted

  if (level=="AG"){
    river$AG[["width"]] <- width_fitted
    river$AG[["depth"]] <- depth_fitted
    river$AG[["discharge"]] <- Q_fitted
    river$AG[["velocity"]] <- velocity_fitted
    river$AG[["volume"]] <- volume
    river$AG[["hydraulicRadius"]] <- RH
    river$AG[["shearStress"]] <- tau
  } else if (level=="RN"){
    river$RN[["width"]] <- width_fitted
    river$RN[["depth"]] <- depth_fitted
    river$RN[["discharge"]] <- Q_fitted
    river$RN[["velocity"]] <- velocity_fitted
    river$RN[["volume"]] <- volume
    river$RN[["hydraulicRadius"]] <- RH
    river$RN[["shearStress"]] <- tau
  }

  invisible(river)
}

Q_GS_generic <- function(depth,ks,slope,width,alpha){
  RH <- eval_RH(width, depth, alpha)
  Q <- ks*width/(1+alpha)*depth*RH^(2/3)*(slope)^0.5
  invisible(Q)
}

depth_GS_generic <- function(Q,ks,slope,width,alpha){
  depth <- (Q/ks/slope^0.5/width)^(3/5) # initial guess with RH=depth
  depth_new <- numeric(length(depth))
  k <- 0
  while (max(abs(depth - depth_new)/depth) > 1e-3){
    k <- k+1
    if (k>1){depth <- depth_new}
    RH <- eval_RH(width, depth, alpha)
    depth_new <- Q/ks/slope^0.5/width/RH^(2/3)*(1+alpha)
  }
  invisible(depth)
}

eval_RH <- function(width, depth, alpha){
  A <- width/(1+alpha)*depth
  w0 <- width/2*depth^(-alpha)
  P <- numeric(length(width))
  for (i in 1:length(P)){
    I <- integrate(function(x) {sqrt(1+ alpha^2*w0[i]^2*x^(2*(alpha-1)))}, 0, depth[i]) # half wetted perimeter
    P[i] <- 2*I$value
    if (P[i]==2*depth[i]){P[i] <- P[i] + width[i]} # correction for rectangular section
  }
  RH <- A/P
  invisible(RH)
}
