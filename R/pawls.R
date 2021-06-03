#' Penalized adaptive weighted least squares regression
#' 
#' Compute weighted least squares regression with \eqn{L_{1}}{L1} regularization on both the
#' coefficients and weight vectors.
#' 
#' @param x a numeric matrix containing the predictor variables without an intercept.  \code{pawls}
#' standardizes the data and includes an intercept by default.
#' @param y a numeric vector containing the response variable.
#' @param nlambda1 the number of lambda1 values (the default is 100).
#' @param nlambda2 the number of lambda2 values (the default is 50).
#' @param lambda1  a numeric vector of non-negative values to be used as tuning parameters of the penalty
#'  for coefficients. By default, a sequence of values of length \code{nlambda1} is computed, equally
#' spaced on the log scale. 
#' @param lambda2  a numeric vector of non-negative values to be used as tunning parameters of the penalty
#'  for weight vectors. By default, a sequence of values of length \code{nlambda2} is computed, equally
#' spaced on the log scale. 
#' @param lambda1.min a numeric value giving the smallest value for \code{lambda1}, as a fraction of the
#' \code{lambda1} maximum. The default is .05. Note that the \code{lambda1}  maximum  is an estimate of tuning parameter that set all the coefficientes to 0.
#' @param lambda2.min a numeric value giving the smallest value for \code{lambda2}, as a fraction of the
#' \code{lambda2} maximum.  The default is .001. Note that the \code{lambda2}  maximum  is an estimate of tuning parameter that set all the weights to 1.
#' @param beta0 the initial estimates of coefficients to be used in the adaptive penalty for \code{beta}.
#' @param w0 the initial estimates of weight vector to be used in the adaptive penalty for  \code{w}.
#' @param initial a character string specifying the initial estimates of both coeffcients and weight vectors in the adaptive penalties.
#' If "\code{uniform}", a non-adaptive \code{pawls} is performed. If "\code{pawls}", then the estimates are obtained by non-adaptive pawls (the 
#' default is "\code{uniform}").
#' @param delta a small positive numeric value as a convergence threshold. The algorithm iterates until the RMSD for the change in both coefficents
#' and weight vectors is less than delta  (the default is 1e-06).
#' @param maxIter a positive numeric value used to determin the maximum number of iterations (the default is 1000).
#' @param intercept  a logical indicating whether a constant term should be 
#' included in the model (the default is \code{TRUE}).
#' @param standardize a logical indicating whether the predictor variables should be normalized 
#' to have unit L2 norm (the default is TRUE). 
#' @param search a character string specifying the algorithm to select tunning parameters for both coefficients and weight
#' vectors. If "cross", the optimal tuning parameters are searched alternatively by minimizing \code{BIC}. If "grid", 
#' the optimal tuning parameters are selected as the pair that minimizs \code{BIC} over a fine grid. The default is "grid".
#' @return 
#' An object of class "pawls.cross"(\code{search=cross}) or \code{"pawls.grid"} 
#' (\code{search=grid}) contaning:
#' 
#' \item{beta}{ a numeric vector containing the respective coefficient estimates with the optimal tuning parameters.}
#' \item{w}{ a numeric vector containing the respective weight estimates with the optimal tuning parameters.}
#' \item{lambda1}{same as above.}
#' \item{lambda2}{same as above.}
#' \item{opt.lambda1}{a numeric value giving the optimal \code{lambda1} in the sense of minimzing \code{BIC}.}
#' \item{opt.lambda2}{a numeric value giving the optimal \code{lambda2} in the sense of minimzing \code{BIC}.}
#' \item{iter}{a numeric matrix with \code{nlambda2} rows and \code{nlambda1} columns giving the number of iteration until convergence
#' at each pair of tuning parameters(\code{search=grid}) or a numeric value giving the number of iteration 
#' in corss search until convergence(\code{search=cross}).}
#' \item{betas}{ a 3-dimension numeric array containing the coefficient estimates. The dimensions are equal to \code{nlambda2}, \code{nlambda1} 
#' and the number of coefficients, respectively. It belongs to the object of class \code{"pawls.grid"} only.}
#' \item{ws}{ a 3-dimension numeric array containing the weight estimates. The dimensions are equal to \code{nlambda2}, \code{nlambda1} 
#' and the number of observations, respectively. It belongs to the object of class \code{"pawls.grid"} only.}
#' \item{raw.bic}{a numeric matrix with \code{nlambda2} rows and \code{nlambda1} columns giving the raw \code{BIC}
#' value for each pair of tuning prameters. It belongs to the object of class \code{"pawls.grid"} only.}
#' \item{bic}{a numeric matrix with \code{nlambda2} rows and \code{nlambda1} columns giving the \code{BIC}
#' value for each pair of tuning prameters. It belongs to the object of class \code{"pawls.grid"} only.}
#' @author Bin Luo, Xiaoli Gao
#' @examples 
#' ## generate data
#' library("mvtnorm")
#' set.seed(123)  # for reproducibility
#' n = 100 # number of observations
#' p = 8 # number of variables
#' beta = c(1, 2, 3, 0, 0, 0, 0, 0) # coefficients
#' sigma <- 0.5      # controls signal-to-noise ratio
#' epsilon <- 0.1    # contamination level
#' Sigma <- 0.5^t(sapply(1:p, function(i, j) abs(i-j), 1:p))
#' x <- rmvnorm(n, sigma=Sigma)    # predictor matrix
#' e <- rnorm(n)                   # error terms
#' i <- 1:ceiling(epsilon*n)       # observations to be contaminated
#' e[i] <- e[i] + 5                # vertical outliers
#' y <- c(x %*% beta + sigma * e)  # response
#' x[i,] <- x[i,] + 5              # bad leverage points
#' 
#' ## fit pawls model over a find grid of tuning parameters
#' pawls(x,y)
#' 
#' ## fit adaptive pawls model over a find grid of tuning parameters
#' pawls(x,y,lambda1.min = 0.001, lambda2.min = 0.05, initial = "PAWLS")
#' 
#' ## fit adaptive pawls model using corss search over a grid of tuning parameters
#' pawls(x,y,lambda1.min = 0.001, lambda2.min = 0.05, initial = "PAWLS", search="cross")
#' @export
pawls = function(x, y, nlambda1 = 100, nlambda2 = 50, lambda1 = NULL, lambda2 = NULL, lambda1.min=0.001,
    lambda2.min=0.01, penalty.factor1 = rep(1,p), penalty.factor2 = rep(1,n), delta = 1e-06, 
    maxIter = 1e03, intercept = TRUE) {
    n = length(y)
    p = dim(x)[2]
    ## check error
    if (!"matrix" %in% class(x)) {
        tmp <- try(x <- as.matrix(x), silent = TRUE)
        if (class(tmp)[1] == "try-error") 
            stop("x must be a matrix or able to be coerced to a matrix")
    }
    if (!"numeric" %in%class(y)) {
        tmp <- try(y <- as.numeric(y), silent = TRUE)
        if (class(tmp)[1] == "try-error") 
            stop("y must numeric or able to be coerced to numeric")
    }
    if (any(is.na(y)) | any(is.na(x))) 
    stop("Missing data (NA's) detected.Take actions to eliminate missing data before passing 
         X and y to pawls.")
    
    ## default setting
    penalty2 <- "L1"
    penalty1 <- "L1"
    
    if (!is.null(lambda1)) 
      nlambda1 <- length(lambda1)
    if (!is.null(lambda2)) 
        nlambda2 <- length(lambda2)


    ## check intercept
    if (intercept) {
      x <- cbind(rep(1, n), x)
      p <- p+1
      peantly.factor1 <- c(0, penalty.factor1)
    }
    
    ## set tunning parameter
    if (is.null(lambda1)||is.null(lambda2)) {
        lambda = setup_parameter(x=x, y=y, nlambda1=nlambda1, nlambda2=nlambda2, 
                                 lambda1.min=lambda1.min, lambda2.min=lambda2.min,
                                 penalty.factor1, penalty.factor2)
        if (is.null(lambda1)) 
          lambda1 = lambda$lambda1
        if (is.null(lambda2)) 
            lambda2 = lambda$lambda2
    }
    
    ## Fit
    
    startBeta = rep(0, p)
    startW = rep(1, n)
    res <- .Call("PAWLS_GRID", x, y, penalty1, penalty2, lambda1, lambda2, penalty.factor1, penalty.factor2, delta, maxIter, 
                 ifelse(intercept, 1, 0), startBeta = startBeta, startW = startW)
    
    fit = list(beta = array(res[[1]], dim = c(nlambda2, nlambda1, p)), w = array(res[[2]], dim = c(nlambda2, nlambda1, n)), 
               wloss = array(res[[3]], dim = c(nlambda2, nlambda1)), loss = array(res[[4]], dim = c(nlambda2, nlambda1)), 
               iter = array(res[[5]], dim = c(nlambda2, nlambda1)))
    fit$lambda1 <- lambda1
    fit$lambda2 <- lambda2
    fit$intercept <- intercept
    fit
}


