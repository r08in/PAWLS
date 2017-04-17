
pawls2 = function(x, y, lambda1 = NULL, lambda2 = NULL, nlambda1 = 50, nlambda2 = 100, 
                 lambda1.min=1e-03, lambda2.min=0.05, beta0 = NULL, w0 = NULL,initial = c("uniform","PAWLS"), 
                 delta = 1e-06, maxIter = 1000, intercept = TRUE, standardize = TRUE, search = c("cross", "all")) {
  
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
  
  penalty1 <- "1-w0"
  penalty2 <- "LASSO"
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
    init = pawls(x, y, intercept = intercept,search = "all")
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
  if (missing(lambda1) || missing(lambda2)||is.null(lambda1)||is.null(lambda2)) {
    lambda = setup_parameter(XX, yy, nlambda1, nlambda2, lambda1.min=lambda1.min, lambda2.min=lambda2.min, beta0, w0, intercept = intercept, penalty1 = penalty1)
    if (is.null(lambda1)) 
      lambda1 = lambda$lambda1
    if (is.null(lambda2)) 
      lambda2 = lambda$lambda2
  }
  
  ## Fit
  if (search == "all") { # search for the whole grid
    res = pawls_grid(XX, yy, penalty1 = penalty1, penalty2 = penalty2, lambda1, lambda2, beta0, w0, delta, 
                     maxIter, intercept = intercept)
    res = BIC_grid(res$wloss, res$beta, res$w, lambda1, lambda2, criterion = criterion)
  } else {# serach="cross"
    res = pawls_cross(XX, yy, penalty1 = penalty1, penalty2 = penalty2, lambda1, lambda2, beta0, w0, delta, maxIter, 
                      intercept = intercept, criterion = criterion, startBeta = startBeta, startW = startW)
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
