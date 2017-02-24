
srcdreg = function(x, y, penalty1 = c("1-w0", "log", "null"), penalty2 = c("LASSO", "RIDGE", "MCP"), lambda1 = NULL, 
    lambda2 = NULL, nlambda1 = 50, nlambda2 = 100, lambda1.min=1e-03, 
    #lambda2.min=0.05, 
    lambda2.min=0.05,
    beta0 = NULL, w0 = NULL, startBeta = NULL, startW = NULL, initial = c("uniform", 
        "LTS", "LASSO", "PAWLS"), delta = 1e-06, maxIter = 1000, intercept = TRUE, standardize = TRUE, updateInitialTimes = 0, 
    criterion = c("BIC", "AIC", "CV"), initCrit=c("BIC", "AIC", "CV"), search = c("cross", "all", "fixw", "crossDynamic"), ...) {
    ## error checking
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
    
    penalty1 <- match.arg(penalty1)
    penalty2 <- match.arg(penalty2)
    initial <- match.arg(initial)
    criterion <- match.arg(criterion)
    initCrit <- match.arg(initCrit)
    search <- match.arg(search)
    
    
    # if (nlambda1 < 2||nlambda2<2) stop('nlambda must be at least 2')
    if (!is.null(lambda1)) 
        nlambda1 = length(lambda1)
    if (!is.null(lambda2)) 
        nlambda2 = length(lambda2)
    
    if (any(is.na(y)) | any(is.na(x))) 
        stop("Missing data (NA's) detected.Take actions to eliminate missing data before passing X and y to srcdreg.")
    
    # initial
    n = length(y)
    p = dim(x)[2]
    if (initial == "LTS") {
        require(robustHD)
        init = sparseLTS(x, y, intercept = intercept)
        beta0 = SetBeta0(init$coefficients)
        w0 = ifelse(init$wt == 1, 0.99, 0.01)
        # w0=UpdateWeight(x,y,beta0) w0=ifelse(as.vector(w0)==1,0.99,w0)
    } else if (initial == "LASSO") {
        init = srcdreg(x, y, penalty1 = "null", nlambda1 = 1, intercept = intercept, search = "fixw")
        beta0 = SetBeta0(init$beta)
        # if(intercept) beta0=c(1,beta0)
        w0 = rep(0.99, n)
    } else if (initial == "PAWLS") {
        init = srcdreg(x, y, intercept = intercept,criterion = initCrit)
        beta0 = SetBeta0(init$beta)
        w0 = ifelse(init$w == 1, 0.99, init$w)
        
        # seting starting beta and w0 require(robustHD) #init_LTS=sparseLTS(x,y,intercept=intercept)
        # load('data\\init_LTS.rda') startBeta=init_LTS$coefficients startW=UpdateWeight(x,y,init_LTS$coefficients)
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
    
    
    # intercept
    if (intercept) {
        x = AddIntercept(x)
    }
    
    # sandardize
    std = 0
    scale = 0
    if (standardize) {
        std <- .Call("GroupStandardize", x, y)
        XX <- std[[1]]
        yy <- std[[2]]
        scale <- std[[3]]
    } else {
        XX = x
        yy = y
    }
    
    ## setup parameter
    if (missing(lambda1) || missing(lambda2)||is.null(lambda1)||is.null(lambda2)) {
        lambda = SetupParameter(XX, yy, nlambda1, nlambda2, lambda1.min=lambda1.min, lambda2.min=lambda2.min, beta0, w0, intercept = intercept, penalty1 = penalty1)
        if (is.null(lambda1)) 
            lambda1 = lambda$lambda1
        if (is.null(lambda2)) 
            lambda2 = lambda$lambda2
    }
    
    ## Fit
    
    if (criterion == "BIC" || criterion == "AIC") {
        # search for the whole grid
        if (search == "all") {
            res = InnerReg(XX, yy, penalty1 = penalty1, penalty2 = penalty2, lambda1, lambda2, beta0, w0, delta, 
                maxIter, intercept = intercept)
            res = BICPWLQ(res$wloss, res$beta, res$w, lambda1, lambda2, criterion = criterion)
        } else if (search == "fixw") {
            # ALasso
            res = InnerReg(XX, yy, penalty1 = "null", penalty2 = penalty2, lambda1 = 1, lambda2, beta0, w0, delta, 
                maxIter, intercept = intercept)
            # res=RCDReg2(XX, yy,penalty1=penalty1,penalty2='ADL',lambda1=1,lambda2,beta0,w0,delta,
            # maxIter,intercept=intercept,fixW=TRUE)
            res = BICPWLQ(res$wloss, res$beta, res$w, lambda1 = 1, lambda2, alpha = 1, criterion = criterion)
        } else if (search == "crossDynamic") {
            while (TRUE) {
                res = RCDReg3(XX, yy, penalty1 = penalty1, penalty2 = penalty2, lambda1, lambda2, beta0, w0, delta, 
                  maxIter, intercept = intercept, updateInitialTimes = updateInitialTimes, criterion = criterion)
                if ((res$index1 < nlambda1 - 1) && (res$index2 < nlambda2 - 1)) 
                  break
                # lambda2 for beta
                if (res$index1 >= nlambda1 - 1) {
                  lambda1 = c(lambda1[1:(nlambda1 - 2)], logSeq(lambda1[nlambda1 - 1], lambda1[nlambda1], nlambda1))
                  nlambda1 = length(lambda1)
                }
                
                if (res$index2 >= nlambda2 - 10) {
                  # lambda2 for beta
                  lambda2 = c(lambda2[1:(nlambda2 - 2)], logSeq(lambda2[nlambda2 - 1], lambda2[nlambda2], nlambda2))
                  nlambda2 = length(lambda2)
                }
                beta0 = res$beta0
                w0 = res$w0
                
            }
            
        } else {
            # (search=='cross')
            res = RCDReg3(XX, yy, penalty1 = penalty1, penalty2 = penalty2, lambda1, lambda2, beta0, w0, delta, maxIter, 
                intercept = intercept, updateInitialTimes = updateInitialTimes, criterion = criterion, startBeta = startBeta, 
                startW = startW)
        }
        
    } else if (criterion == "CV") {
        res = RCDReg4(XX, yy, penalty1 = penalty1, penalty2 = penalty2, lambda1, lambda2, beta0, w0, delta, maxIter, 
            intercept = intercept)
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
