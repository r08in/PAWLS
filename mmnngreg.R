#MM-non negative garrote regression
mmnngreg=function(x,y)
{
  require(robustbase)
  #Compute initial estimate with MM-regression:
  #beta_init <- coef(lmrob(y~x))
  #beta_init
  
  require(robustHD)
  beta_init=sparseLTS(x,y)$coefficients
  #Compute S-nonnegative garrote estimate using the initial estimate beta_init
  #and the BIC-criterion to select lambda:
  sol_S <- SnngarroteGrid(x,y,TRUE, lambda=seq(0.001,0.5,length=50), crit="BIC", beta=beta_init)
  sol_S
  
  #Compute MM-nonnegative garrote estimate using the initial estimate beta_init,
  #the shrinkage coefficients of sol_S as starting point of the iterative procedure 
  #and the BIC-criterion to select lambda:
  sol_MM <- MMnngarroteGrid(x,y, TRUE, lambdaMM=seq(0.001,0.2,length=50), critMM="BIC", beta=beta_init, snng=list(cinit=sol_S$C, scaleinit=sol_S$scale))
  sol_MM
}

# S-nonnegative garrote method

require(robustbase)
require(Rcpp)
require(parallel)
#search("SNNG.dll")
dyn.load('D:\\RProject\\RobustCD\\SNNG.dll')


# Estimate the coefficients of multiple linear regression model with S-nonnegative garrote method,
# over a grid of values for the penalty parameter.
# 
# 
# X:  a numeric matrix containing the predictor variables.
# Y:  a numeric vector containing the response variable.
# intercept:  a logical indicating whether a constant term should be 
#             included in the model (the default is TRUE).
# lambda:  a numeric vector of non-negative numeric values to be used as 
#          penalty parameter.
# crit:  a character string specifying the optimality criterion to be 
#        used for selecting the final model.  Possible values are "L" for 
#        the L-curve criterion, "BIC" for the robust Bayes information criterion,  
#        and "AIC" for the robust Aikake's information criterion.
# beta:  a numeric vector containing initial estimators for coefficients.
# N:  a numeric value giving the number of initial subsamples to be used in 
#     the first phase of the algorithm. The default is 500 initial subsamples. 
# k:  a positive integer giving the number of I-steps to perform on 
#     all subsamples in the first phase of the algorithm (the default is to 
#     perform two I-steps).
# bestr:  a numeric value giving the number of subsamples to keep 
#         after the first phase of k I-steps. For those remaining subsets, 
#         additional I-steps are performed until convergence. The default is 5 subsamples. 
# b: a numeric value used to obtain consistency for the scale at the normal 
#    error model. The default is 0.5.
# cc: a numeric value used in the loss function of Tukey's biweight family.
#     The default is 1.5476.
# ncores:  a positive integer giving the number of processor cores to be 
#          used for parallel computing (the default is 1 for no parallelization).  If 
#          this is set to NA, all available processor cores are used.
# seed:  optional initial seed for the random number generator.
#
#
#
# Output variables:
#
# betac: a numeric vector of S-nonnegative garrote estimates. 
# scale: the value of the M-scale.
# C:  a numeric vector of the shrinkage factors.
# beta:  a numeric vector of initial estimates.
# lambda: the optimal value for the penalty parameter.


SnngarroteGrid <- function(X, Y, intercept=TRUE, lambda, crit=c("L", "BIC", "AIC"), beta, N = 500, k=2, bestr = 5, b = 0.5, cc = 1.5476, ncores=1, seed = NULL){
  n <- length(Y)
  d <- dim(X)
  if(!isTRUE(n==d[1])) stop(sprintf("'X' must have %d rows", n))
  if(missing(lambda)){
    stop("missing 'lambda'")
  }
  if(!is.numeric(lambda) || length(lambda)==0 || any(!is.finite(lambda))){
    stop("missing or invalid value of 'lambda'")
  }
  if(any(negative <- lambda < 0)){
    lambda[negative] <- 0
    warning("negative value for 'lambda', using no penalization")
  }
  lambda <- sort(unique(lambda), decreasing = TRUE)
  intercept <- isTRUE(intercept)
  if(intercept){
    X <- cbind(rep(1,n), X)
    d[2] <- d[2] + 1
  }
  initial <- replicate(N, sample(n, d[2],replace=TRUE))
  if(missing(beta) || any(!is.finite(beta))){
    beta <- as.vector(coef(lmrob(Y~X-1)))
    warning("missing or infinite values of 'beta'; using default values")
  }
  crit <- match.arg(crit)
  ncores <- rep(ncores, length.out = 1)
  if(is.na(ncores)) ncores <- detectCores()
  if(!is.numeric(ncores) || is.infinite(ncores) || ncores < 1){
    ncores <- 1
    warning("infinite value of 'ncores'; using default value") 
  } else {
    ncores <- as.integer(ncores)
  }
  if(!is.null(seed)) set.seed(seed)
  fit <- lapply(lambda, function(l, X, Y, beta, N, k, bestr, b, cc, ncores, seed){Snngarrote(X, Y, intercept=FALSE, lambda=l, beta,  N, k, bestr, b, cc, ncores, seed)}, X, Y, beta, N, k, bestr, b, cc, ncores, seed)
  if(crit == "BIC"){
    critValues <- sapply(fit, AIC, n=n, k=log(n))
    sOpt <- which.min(critValues)
  } 
  else if(crit == "AIC"){
    critValues <- sapply(fit, AIC, n=n)   
    sOpt <- which.min(critValues)
  }
  else if(crit == "L"){
    critValues <- sapply(fit, Lcurve)
    m1 <- min(critValues[1,])
    M1 <- max(critValues[1,])
    m2 <- min(critValues[2,])
    M2 <- max(critValues[2,])
    critValues[1,] <- (critValues[1,]-m1)/(M1-m1) 
    critValues[2,] <- (critValues[2,]-m2)/(M2-m2)
    sOpt <- corner(critValues[ 1,], critValues[2,]) 
  }
  fit <- list(betac=fit[[sOpt]]$betac, scale=fit[[sOpt]]$scale, C=fit[[sOpt]]$C,beta=beta,  lambda=lambda[sOpt])
  return(fit)
}



# Estimate the coefficients of multiple linear regression model with S-nonnegative garrote method
#
#
# Input variables:
#
# X:  a numeric matrix containing the predictor variables.
# Y:  a numeric vector containing the response variable.
# intercept:  a logical indicating whether a constant term should be 
#             included in the model (the default is TRUE).
# lambda:  a non-negative numeric value giving the penalty parameter.
# beta:  a numeric vector containing initial estimators for coefficients.
# N:  a numeric value giving the number of initial subsamples to be used in 
#     the first phase of the algorithm. The default is 500 initial subsamples. 
# k:  a positive integer giving the number of I-steps to perform on 
#     all subsamples in the first phase of the algorithm (the default is to 
#     perform two I-steps).
# bestr:  a numeric value giving the number of subsamples to keep 
#         after the first phase of k I-steps. For those remaining subsets, 
#         additional I-steps are performed until convergence. The default is 5 subsamples. 
# b: a numeric value used to obtain consistency for the scale at the normal 
#    error model. The default is 0.5.
# cc: a numeric value used in the loss function of Tukey's biweight family.
#     The default is 1.5476.
# ncores:  a positive integer giving the number of processor cores to be 
#          used for parallel computing (the default is 1 for no parallelization).  If 
#          this is set to NA, all available processor cores are used.
# seed:  optional initial seed for the random number generator.
#
#
#
# Output variables:
#
# betac: a numeric vector of S-nonnegative garrote estimates. 
# scale: the value of the M-scale.
# C:  a numeric vector of the shrinkage factors.
# beta:  a numeric vector of initial estimates.

Snngarrote <- function(X, Y, intercept=TRUE, lambda, beta, N = 500, k=2, bestr = 5, b = 0.5, cc = 1.5476, ncores=1, seed = NULL){
  n <- length(Y)
  X <- as.matrix(X)
  d <- dim(X)
  if(!isTRUE(n == d[1])) stop(sprintf("'X' must have %d rows", n))
  if(missing(lambda)){ 
    stop("missing 'lambda'")
  } 
  intercept <- isTRUE(intercept)
  if(intercept){
    X <- cbind(rep(1,n), X)
    d[2] <- d[2] + 1
  }
  if(missing(beta) || any(!is.finite(beta))){
    beta <- as.vector(coef(lmrob(Y~X-1)))
    warning("missing or infinite values of 'beta'; using default values")
  }
  Z <- as.matrix(X%*%diag(beta))
  ncores <- rep(ncores, length.out = 1)
  if(is.na(ncores)) ncores <- detectCores()
  if(!is.numeric(ncores) || is.infinite(ncores) || ncores < 1){
    ncores <- 1
    warning("infinite value of 'ncores'; using default value") 
  } else {
    ncores <- as.integer(ncores)
  }
  if(!is.null(seed)) set.seed(seed)
  if(!is.null(seed)) set.seed(seed)
  initial <- replicate(N, sample(n, d[2],replace=TRUE))
  fit <- .Call("R_Snngarrote", RZ=Z, R_Y=Y, R_lambda=lambda, R_initial=initial, R_k=k, R_bestr=bestr, R_b=b, R_cc=cc, R_ncores=ncores)  
  fit$coefficients <- drop(fit$coefficients)
  betac <- beta * fit$coefficients
  fit <- list(betac= betac, scale = fit$scale,  C = fit$coefficients, beta = beta)
  return(fit)
}




AIC <- function(object, n, k = 2) {
  scale <- object$scale
  df <- length(which(object$betac!=0))
  log(scale) + k*df/n
}

Lcurve <- function(object){
  L <- rep(0,2)
  L[1] = sum(object$C)
  L[2] = object$scale
  return(L)
}


corner <- function(a, b){
  distance <- sqrt(a^2 + b^2)
  which.min(distance)
}

# MM-nonnegative garrote method

require(robustbase)
require(Rcpp)
require(parallel)
dyn.load('D:\\RProject\\RobustCD\\SNNG.dll')


# Estimate the coefficients of multiple linear regression model with MM-nonnegative garrote method,
# over a grid of values for the penalty parameter.
# 
# 
# X:  a numeric matrix containing the predictor variables.
# Y:  a numeric vector containing the response variable.
# intercept:  a logical indicating whether a constant term should be 
#             included in the model (the default is TRUE).
# lambdaMM:  a numeric vector of non-negative numeric values to be used as 
#            penalty parameter.
# critMM:  a character string specifying the optimality criterion to be 
#        used for selecting the final model.  Possible values are "L" for 
#        the L-curve criterion, "BIC" for the robust Bayes information criterion,  
#        and "AIC" for the robust Aikake's information criterion.
# beta:  a numeric vector containing initial estimators for coefficients.
# ccMM: a numeric value used in the loss function of Tukey's biweight family.
#       The default is 5.182.
# snng: list with either cinit and scaleinit or lambdaS and critS.
#       -cinit: a numeric vector containing initial estimates for shrinkage factors.
#       -scaleinit: a numeric value giving the estimate for scale.
#       -critS: if missing cinit, selection criterium for lambdaS.
#       -lambdaS:  if missing cinit, vector of regularization parameters for S-NNG.
#
#
#
# Output variables:
#
# betac: a numeric vector of MM-nonnegative garrote estimates. 
# scale: the value of the M-scale.
# C:  a numeric vector of the shrinkage factors.
# beta:  a numeric vector of initial estimates.
# lambda: the optimal value for the penalty parameter.


MMnngarroteGrid <- function(X, Y, intercept, lambdaMM, critMM=c("L", "BIC", "AIC"),  beta, ccMM=5.182, snng=list(cinit, scaleinit, lambdaS, critS)){
  n <- length(Y)
  d <- dim(X)
  if(!isTRUE(n==d[1])) stop(sprintf("'X' must have %d rows", n))
  if(missing(lambdaMM)){
    stop("missing 'lambdaMM'")
  }
  lambdaMM <- sort(unique(lambdaMM), decreasing = TRUE)
  intercept <- isTRUE(intercept)
  if(intercept){
    X <- cbind(rep(1,n), X)
    d[2] <- d[2] + 1
  }
  if(missing(beta) || any(!is.finite(beta))){
    beta <- as.vector(coef(lmrob(Y~X-1)))
  }
  critMM <- match.arg(critMM)
  ccMM <- rep(ccMM, length.out = 1)
  if(!is.numeric(ccMM) || !is.finite(ccMM)){
    ccMM <- 5.182
  }
  if(missing(snng)){
    stop("no initial values for shrinkage factors and scale estimate")
  }
  if(is.null(snng$cinit) || any(!is.finite(snng$cinit))){
    temp <- SnngarroteGrid(X,Y,intercept=FALSE, lambda=snng$lambdaS, crit=snng$critS, beta=beta)
    snng$cinit=as.vector(temp$C)
    snng$scaleinit = temp$scale
  }
  Z <- as.matrix(X%*%diag(beta))
  snng$scaleinit <- rep(snng$scaleinit, length.out = 1)
  if(!is.numeric(snng$scaleinit) || !is.finite(snng$scaleinit)){
    snng$scaleinit = Mscale(Y-Z%*%snng$cinit, 0.5, 1.5476)
  }
  fit <- lapply(lambdaMM, function(l, X, Y, beta, ccMM, snng ){MMnngarrote(X, Y,  intercept=FALSE, lambdaMM=l, beta, ccMM, snng)}, X, Y, beta, ccMM, snng)
  if(critMM == "BIC"){
    critValues <- sapply(fit, AICMM,X,Y, ccMM, n, log(n))
    sOpt <- which.min(critValues)
  } 
  else if(critMM == "AIC"){
    critValues <- sapply(fit, AICMM, X, Y, ccMM, n)   
    sOpt <- unname(which.min(critValues))
  }
  else if(critMM == "L"){
    critValues <- sapply(fit, LcurveMM,X,Y, cc=ccMM)
    m1 <- min(critValues[1,])
    M1 <- max(critValues[1,])
    m2 <- min(critValues[2,])
    M2 <- max(critValues[2,])
    critValues[1,] <- (critValues[1,]-m1)/(M1-m1) 
    critValues[2,] <- (critValues[2,]-m2)/(M2-m2)
    sOpt <- corner(critValues[ 1,], critValues[ 2,]) 
  }
  fit <- list(betac=fit[[sOpt]]$betac, scale=fit[[sOpt]]$scale, C=fit[[sOpt]]$C,beta=beta,  lambda=lambdaMM[sOpt])
  fit
}


# Estimate the coefficients of multiple linear regression model with MM-nonnegative garrote method
#
#
# Input variables:
#
# X:  a numeric matrix containing the predictor variables.
# Y:  a numeric vector containing the response variable.
# intercept:  a logical indicating whether a constant term should be 
#             included in the model (the default is TRUE).
# lambdaMM:  a non-negative numeric value giving the penalty parameter.
# beta:  a numeric vector containing initial estimators for coefficients.
# ccMM: a numeric value used in the loss function of Tukey's biweight family.
#       The default is 5.182.
# snng: list with either cinit and scaleinit or lambdaS and critS.
#       -cinit: a numeric vector containing initial estimates for shrinkage factors.
#       -scaleinit: a numeric value giving the estimate for scale.
#       -critS: if missing cinit, selection criterium for lambdaS.
#       -lambdaS:  if missing cinit, vector of regularization parameters for S-NNG.
#
#
#
# Output variables:
#
# betac: a numeric vector of MM-nonnegative garrote estimates. 
# scale: the value of the M-scale.
# C:  a numeric vector of the shrinkage factors.
# beta:  a numeric vector of initial estimates.


MMnngarrote <- function(X, Y, intercept=TRUE, lambdaMM, beta,  ccMM = 5.182, snng=list(cinit, scaleinit, lambdaS, critS)){
  n <- length(Y)
  X <- as.matrix(X)
  d <- dim(X)
  if(!isTRUE(n == d[1])) stop(sprintf("'X' must have %d rows", n))
  lambdaMM <- rep(lambdaMM, length.out = 1)
  if(!is.numeric(lambdaMM) || !is.finite(lambdaMM)) {
    stop("missing or invalid value of 'lambdaMM'")
  }
  intercept <- isTRUE(intercept)
  if(intercept){
    X <- cbind(rep(1,n), X)
    d[2] <- d[2] + 1
  }
  if(missing(beta) || any(!is.finite(beta))){
    beta <- as.vector(coef(lmrob(Y~X-1)))
    warning("missing or infinite values of 'beta'; using default values")
  }
  Z <- as.matrix(X%*%diag(beta))
  if(missing(snng)){
    stop("no initial values for shrinkage factors and scale estimate")
  }
  if(is.null(snng$cinit) || any(!is.finite(snng$cinit))){
    temp <- SnngarroteGrid(X,Y,intercept=FALSE, lambda=snng$lambdaS, crit=snng$critS, beta=beta)
    snng$cinit=as.vector(temp$C)
    snng$scaleinit = temp$scale
  }
  snng$scaleinit <- rep(snng$scaleinit, length.out = 1)
  if(!is.numeric(snng$scaleinit) || !is.finite(snng$scaleinit)){
    snng$scaleinit = Mscale(Y-Z%*%snng$cinit, 0.5, 1.5476)
  }
  ccMM <- rep(ccMM, length.out = 1)
  if(!is.numeric(ccMM) || !is.finite(ccMM)){
    ccMM <- 5.182
    warning("missing of infinite value of 'ccMM'; using default value") 
  }
  fit <- .Call("R_MMnngarrote", RZ=Z, R_Y=Y, R_lambda=lambdaMM, R_initcoef=snng$cinit, R_initscale=snng$scaleinit,R_cc=ccMM)
  betac <- beta * fit$coefficients
  fit <- list(betac= betac, scale = fit$scale, C = fit$coefficients, beta = beta)
  return(fit)
}


AICMM <- function(object,X,Y, cc,n, k = 2) {
  scale <- object$scale
  residuals <- Y-X%*%object$betac
  df <- length(which(object$betac!=0))
  # compute AIC 
  log(mean(rho(residuals/scale, cc))) + k*df/n
}

LcurveMM <- function(object, X,Y, cc){
  L <- rep(0,2)
  residuals <- Y-X%*%object$betac
  L[1] = sum(object$C)
  L[2] = mean(rho(residuals/object$scale, cc))
  return(L)
}


rho <- function(u, cc){
  w <- abs(u)<=cc
  v <- (u^2/(2)*(1-(u^2/(cc^2))+(u^4/(3*cc^4))))*w +(1-w)*(cc^2/6)
  v <- v*6/cc^2
  return(v)
}

Mscale <- function(u, b, cc, initial.sc=median(abs(u))/.6745) {
  # find the scale, full iterations
  max.it <- 200
  # magic number alert
  #sc <- median(abs(u))/.6745
  sc <- initial.sc
  i <- 0
  eps <- 1e-20
  # magic number alert
  err <- 1
  while( ( (i <- i+1) < max.it ) && (err > eps) ) {
    sc2 <- sqrt( sc^2 * mean( rho( u / sc, cc ) ) / b   )
    err <- abs(sc2/sc - 1)
    sc <- sc2
  }
  return(sc)
}

corner <- function(a, b){
  distance <- sqrt(a^2 + b^2)
  which.min(distance)
}
