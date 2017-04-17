#' Penalized adaptive weighted least squares regression
#' 
#' Compute weighted least squares regression with \eqn{L_{1}}{L1} regularization on both the
#' coefficients and weight vectors.
#' 
#' @param x a numeric matrix containing the predictor variables without an intercept.  \code{pawls}
#' standardizes the data and includes an intercept by default.
#' @param y a numeric vector containing the response variable.
#' @param nlambda1 The number of lambda1 values (the default is 100).
#' @param nlambda2 The number of lambda2 values (the default is 50).
#' @param lambda1  a numeric vector of non-negative values to be used as penalty parameter for coefficients. 
#' By default, a sequence of values of length \code{nlambda1} is computed, equally
#' spaced on the log scale. 
#' @param lambda2  a numeric vector of non-negative values to be used as penalty parameter for weight vectors. 
#' By default, a sequence of values of length \code{nlambda2} is computed, equally
#' spaced on the log scale. 
#' @param lambda1.min a numeric value giving the ratio of minimum \code{lambda1} to maximum \code{lambda1}. 
#' The maximum \code{lambda1} is an estimate of penalty parameter that set all the coefficientes to 0.
#' @param lambda2.min a numeric value giving the ratio of minimum \code{lambda2} to maximum \code{lambda2}. 
#' The maximum \code{lambda2} is an estimate of penalty parameter that set all the weight to 1.
#' @param beta0 the initial estimates of coefficients \code{beta} used in the adaptive penalty.
#' @param w0 the initial estimates of weight vector \code{w} used in the adaptive penalty.
#' @param initial a character string specifying the initial estimates of both coeffcients and weight vectors in the adaptive penalties.
#' If "\code{uniform}", a non-adaptive \code{pawls} is performed. If "\code{pawls}", then the estimates are obtained by non-adaptive pawls (the 
#' default is "\code{uniform}").
#' @param delta a small positive numeric value used to determine whether the variability within a variable is 
#' too small (the default is 1e-06).
#' @param maxIter a positive numeric value used to determin the maximum number of iteration for optimization.
#' @param intercept  a logical indicating whether a constant term should be 
#' included in the model (the default is \code{TRUE}).
#' @param standardize a logical indicating whether the predictor variables should be normalized 
#' to have unit L2 norm (the default is TRUE). 
#' @param search a character string specifying the algorithm to select tunning parameters for both coefficients and weight
#' vectors. If "cross", the optimal tuning parameters are searched alternatively by minimizing \code{BIC}. If "grid", 
#' the optimal tuning parameters are selected as the pair that minimiz \code{BIC} over a fine grid.
pawls = function(x, y, nlambda1 = 100, nlambda2 = 50, lambda1 = NULL, lambda2 = NULL, lambda1.min=0.05,
    lambda2.min=1e-03, beta0 = NULL, w0 = NULL,initial = c("uniform","PAWLS"), delta = 1e-06, 
    maxIter = 1e03, intercept = TRUE, standardize = TRUE, search = c("cross", "grid")) {
  
    ## check error
    if (class(x) != "matrix") {
        tmp <- try(x <- as.matrix(x), silent = TRUE)
        if (class(tmp)[1] == "try-error") 
            stop("x must be a matrix or able to be coerced to a matrix")
    }
    if (class(y) != "numeric") {
        tmp <- try(y <- as.numeric(y), silent = TRUE)
        if (class(tmp)[1] == "try-error") 
            stop("y must numeric or able to be coerced to numeric")
    }
    if (any(is.na(y)) | any(is.na(x))) 
    stop("Missing data (NA's) detected.Take actions to eliminate missing data before passing 
         X and y to pawls.")
    initial <- match.arg(initial)
    search <- match.arg(search)
    
    ## default setting
    penalty2 <- "L1"
    penalty1 <- "L1"
    criterion <- "BIC"
    startBeta <- NULL
    startW <- NULL
    
    if (!is.null(lambda1)) 
      nlambda1 <- length(lambda1)
    if (!is.null(lambda2)) 
        nlambda2 <- length(lambda2)

    ## set initial
    n = length(y)
    p = dim(x)[2]
    if (initial == "PAWLS") {
      init = pawls(x, y, intercept = intercept,search = "grid")
      beta0 = SetBeta0(init$beta)
      w0 = ifelse(init$w == 1, 0.99, init$w)
    } else if (initial == "uniform") {
      if (is.null(beta0)) {
        beta0 = rep(1, p)
        if (intercept) 
          beta0 = c(1, beta0)
      }
      if (is.null(w0)) 
        w0 = rep(0.99, n)
      beta0 = SetBeta0(beta0)
      w0 = ifelse(w0 == 1, 0.99, w0)
    }
    
    ## check intercept
    if (intercept) {
        x = AddIntercept(x)
    }
    
    ## sandardize
    std = 0
    scale = 0
    if (standardize) {
        std <- .Call("Standardize", x, y)
        XX <- std[[1]]
        yy <- std[[2]]
        scale <- std[[3]]
    } else {
        XX = x
        yy = y
    }
    
    ## set tunning parameter
    if (is.null(lambda1)||is.null(lambda2)) {
        lambda = setup_parameter(x=XX, y=yy, nlambda1=nlambda1, nlambda2=nlambda2, 
                                 lambda1.min=lambda1.min, lambda2.min=lambda2.min, beta0=beta0, w0=w0)
        if (is.null(lambda1)) 
          lambda1 = lambda$lambda1
        if (is.null(lambda2)) 
            lambda2 = lambda$lambda2
    }
    
    ## Fit
    if (search == "grid") { # search for the whole grid
      res = pawls_grid(x=XX, y=yy, penalty1 = penalty1, penalty2 = penalty2, lambda1=lambda1, lambda2=lambda2,
                     beta0=beta0, w0=w0, delta=delta, maxIter=maxIter, intercept = intercept)
      res = BIC_grid(res$wloss, res$beta, res$w, lambda2, lambda1, criterion = criterion)
      class(res) <- "grid.pawls"
    } else {# serach="cross"
      res = pawls_cross(x=XX, y=yy, penalty1 = penalty1, penalty2 = penalty2, lambda1=lambda1, lambda2=lambda2, 
                        beta0=beta0, w0=w0, delta=delta, maxIter=maxIter, intercept = intercept, 
                        criterion = criterion, startBeta = startBeta, startW = startW)
      class(res) <- "cross.pawls"
    }
    
    ## unstandardize
    if (standardize) {
        scale = ifelse(scale == 0, 0, 1/scale)
        res$beta = res$beta * scale
    }
    res
}

SetBeta0 = function(beta0) {
    b = 1e-6
    ifelse(abs(beta0) < b, b, beta0)
}

AddIntercept = function(x) {
    n = dim(x)[1]
    x1 = rep(1, n)
    xx = cbind(x1, x)
    xx
}
