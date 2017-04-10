pamls= function(x, y, penalty1 = c("L1"), penalty2 = c("L1"), lambda1 = NULL, 
    lambda2 = NULL, nlambda1 = 50, nlambda2 = 100, 
    lambda1.min=1e-03, 
    lambda2.min=0.05, 
    beta0 = NULL, gam0 = NULL, startBeta = NULL, startGam = NULL, initial = "uniform", delta = 1e-06, maxIter = 1000, intercept = TRUE, standardize = TRUE, updateInitialTimes = 0, 
    criterion = c("BIC", "AIC"), initCrit=c("BIC", "AIC"), search = c("cross", "all"), ...) {
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
    
    criterion <- match.arg(criterion)
    initCrit <- match.arg(initCrit)
    search <- match.arg(search)
    
    
    # if (nlambda1 < 2||nlambda2<2) stop('nlambda must be at least 2')
    if (!is.null(lambda1)) 
        nlambda1 = length(lambda1)
    if (!is.null(lambda2)) 
        nlambda2 = length(lambda2)
    
    if (any(is.na(y)) | any(is.na(x))) 
        stop("Missing data (NA's) detected.Take actions to eliminate missing data before passing X and y to pamls.")
    
    # initial
    n = length(y)
    p = dim(x)[2]
    if (initial == "PAMLS") {
      init = pamls(x, y, intercept = intercept,criterion = initCrit, search = "all")
      beta0 = init$beta
      gam0 = init$gam
    } else if (initial == "uniform") {
        if (is.null(beta0)) {
            beta0 = rep(1, p)
            if (intercept) 
                beta0 = c(1, beta0)
        }
        if (is.null(gam0)) {
            gam0 = rep(1, n)
        }
    }
    b = 1e-6  
    beta0 = ifelse(abs(beta0) < b, b, beta0)
    gam0 = ifelse(abs(gam0) < b, b, gam0)
    
    
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
        lambda = ComputeParameter(XX, yy, nlambda1, nlambda2, lambda1.min=lambda1.min, lambda2.min=lambda2.min, beta0, gam0)
        if (is.null(lambda1)) 
            lambda1 = lambda$lambda1
        if (is.null(lambda2)) 
            lambda2 = lambda$lambda2
    }
    
    ## Fit
    
    if (criterion == "BIC" || criterion == "AIC") {
        # search for the whole grid
        if (search == "all") {
            res = Innerpamls(XX, yy, lambda1, lambda2, beta0, gam0, delta, 
                maxIter, intercept = intercept)
            res = BIC4PAMLS(res$loss, res$beta, res$gam, lambda1, lambda2, criterion = criterion)
        } else {
            # (search=='cross')
            res = pamlsreg(XX, yy, lambda1, lambda2, beta0, gam0, delta, maxIter, 
                intercept = intercept, criterion = criterion, startBeta = startBeta, startGam = startGam)
        }
        
    } 
    ## unstandardize
    if (standardize) {
        scale = ifelse(scale == 0, 0, 1/scale)
        res$beta = res$beta * scale
    }
    res
}
# without search the whole grid of lambda
pamlsreg= function(x, y,lambda1, lambda2, beta0, gam0, delta, maxIter, intercept = TRUE,
                   criterion = "BIC", startBeta = NULL, startGam = NULL) {
  L1 = length(lambda1)
  L2 = length(lambda2)
  m = dim(x)[2]
  n = length(y)
  pre1 = pre2 = 0
  lmaxIter = 50
  
  # initial lambda
  index2 = round(L2/2)
  res = Innerpamls(x, y, lambda1, lambda2[index2], beta0, gam0, delta, maxIter, 
                 intercept = intercept, startBeta = startBeta, startGam = startGam)  #fix lambda2
  bic1 = BIC4PAWLS(as.vector(res$loss), apply(matrix(res$gam, L1, n) != 0 + 0, 1, sum), 
                   apply(matrix(res$beta, L1, m) != 0 + 0, 1, sum), n, m, criterion = criterion) 
  index1 =bic1$index  #find best lambda1
  crit1 = bic1$crit
  
  res = Innerpamls(x, y,lambda1[index1], lambda2, beta0, gam0, delta, maxIter, 
                 intercept = intercept, startBeta = startBeta, startGam = startGam)  #fix lambda1
  bic2 = BIC4PAWLS(as.vector(res$loss), apply(matrix(res$gam, L2, n) != 0 + 0, 1, sum),
                   apply(matrix(res$beta,L2, m) != 0 + 0, 1, sum), n, m, 
                   criterion = criterion)
  index2 = bic2$index  #find best lambda2
  crit2 = bic2$crit
  ## loop to estimate and find the best
  iter = 0
  
  while ((pre1 != index1 || pre2 != index2) && (iter < lmaxIter)) {
    iter = iter + 1
    pre1 = index1
    pre2 = index2
    
    resg = Innerpamls(x, y, lambda1, lambda2[index2], beta0, gam0, delta, 
                    maxIter, intercept = intercept, startBeta = startBeta, startGam = startGam)  #fix lambda2
    bic1 = BIC4PAWLS(as.vector(resg$loss), apply(matrix(resg$gam, L1, n) != 0 + 0, 1, sum), 
                     apply(matrix(resg$beta, L1, m) != 0 + 0, 1, sum), n, m, criterion = criterion) 
    index1 =bic1$index  #find best lambda1
    crit1 = bic1$crit
    
    res = Innerpamls(x, y,lambda1[index1], lambda2, beta0, gam0, delta, 
                   maxIter, intercept = intercept, startBeta = startBeta, startGam = startGam)  #fix lambda1
    bic2 = BIC4PAWLS(as.vector(res$loss), apply(matrix(res$gam, L2, n) != 0 + 0, 1, sum),
                     apply(matrix(res$beta,L2, m) != 0 + 0, 1, sum), n, m, 
                     criterion = criterion)
    index2 = bic2$index  #find best lambda2
    crit2 = bic2$crit
    if (pre2 == index2 && pre1 == index1) 
    {
      break
    }
    
  }
  
  list(lambda1 = lambda1[index1], lambda2 = lambda2[index2], beta = as.vector(res$beta[1, index2, ]), 
       gam = as.vector(res$gam[1,index2, ]), loss = res$loss[1,index2], bdf = sum(res$beta[1, index2, ] != 0), 
       gdf = sum(res$gam[1, index2,] != 0), index1 = index1, index2 = index2, iter = iter, gam0 = gam0, 
       beta0 = beta0, lambda1s = lambda1, lambda2s = lambda2, 
       betas = matrix(res$beta, L2, m), gams = matrix(resg$gam, L1, n), crit1=crit1, crit2=crit2)
  
}

ComputeParameter = function(x, y, nlambda1, nlambda2,lambda1.min=1e-03, lambda2.min=0.05, beta0, gam0) {
  n = length(y)
  p =dim(x)[2]
  lambda1Max = max(abs(y * gam0 / n))
  lambda2Max = max(abs(t(x) %*% y/n) * abs(beta0))  # max |betaj|*|xj'y/n|
  lambda1 = logSeq2(lambda1Max, lambda1Max * lambda1.min, nlambda1)
  lambda2 = logSeq2(lambda2Max, lambda2Max * lambda2.min, nlambda2)
  return(list(lambda1 = lambda1, lambda2 = lambda2))
}

Innerpamls = function(x, y,lambda1, lambda2, beta0, gam0, delta, maxIter, intercept = TRUE, 
                      startBeta = NULL, startGam = NULL) {
  # return (res2=InnerReg2(x, y,penalty1,penalty2,lambda1,lambda2,beta0,w0,delta,
  # maxIter,intercept,startBeta=startBeta,startW=startW))
  L1 = length(lambda1)
  L2 = length(lambda2)
  m = dim(x)[2]
  n = dim(x)[1]
  if (is.null(startBeta)) {
    startBeta = rep(0, m)
  }
  if (is.null(startW)) {
    startGam = rep(0, n)
  }
  res <- .Call("Innerpamls", x, y,lambda1, lambda2, beta0, gam0, delta, maxIter, ifelse(intercept,1, 0), startBeta = startBeta, startGam = startGam)
  
  res = list(beta = array(res[[1]], dim = c(L1, L2, m)), gam = array(res[[2]], dim = c(L1, L2, n)),
             loss = array(res[[3]], dim = c(L1, L2)), iter = array(res[[4]], dim = c(L1, L2)))
  
  res
  
}

Innerpamls2 = function(x, y,lambda1, lambda2, beta0, gam0, delta, maxIter, intercept = TRUE, 
                     startBeta = NULL, startGam = NULL) {
  
  ## declaration
  n = length(y)
  m = dim(x)[2]
  p = m
  L1 = length(lambda1)
  L2 = length(lambda2)
  lstart1 = 1  ##should be 1
  lstart2 = 1
  
  ## reslut to be returned
  beta = array(0, dim = c(L1, L2, m))
  gam = array(0, dim = c(L1, L2, n))
  loss = matrix(ncol = L2, nrow = L1, 0)
  iter = matrix(ncol = L2, nrow = L1, 0)
  
  betaPre = rep(0, m)
  gamPre = rep(0, n)
  if (!is.null(startBeta)) {
    betaPre = startBeta
  }
  if (!is.null(startGam)) {
    gamPre = startGam
  }
  r = y - x %*% betaPre
  ## iteration for each lamda1
  for (l1 in lstart1:L1) {
    ## initial
    
    shift = rep(0, m + n)
    loss[l1, 1] = t(r - gamPre) %*% (r - gamPre)  ##initial loss[l1,1]
    lam1 = (lambda1[l1]/abs(gam0)) * n
    
    ## iteration for each lamda2
    for (l2 in lstart2:L2) {
      lam2 = lambda2[l2]/abs(beta0)
      if (intercept) {
        lam2[1] = 0
      }
      ## iteration for all covariates
      while (iter[l1, l2] < maxIter) {
        iter[l1, l2] = iter[l1, l2] + 1
        
        ## iteration for each beta
        for (j in 1:m) {
          ## (1)calculate zj
          zj = t(x[, j]) %*% (r - gamPre)/n + betaPre[j]
          ## (2)update betaj
          beta[l1, l2, j] = SoftThreshold(zj, lam2[j])

          ## (3)update r
          shift[j] = beta[l1, l2, j] - betaPre[j]
          r = r - x[, j] * shift[j]
        }
        
        ## update gamma
        gam[l1, l2, ]=SoftThreshold(r,lam1)
        shift[(m + 1):(m + n)] = gam[l1, l2, ] - gamPre
        
        ## update betaPre and gamPre for next iteration
        betaPre = beta[l1, l2, ]
        gamPre = gam[l1, l2, ]
        ## Check for convergence
        if (t(shift) %*% shift < delta) {
          break
        }
        
      }  #end for the inner loop
      
      ## compute square of loss
      loss[l1, l2] = t(r - gamPre) %*% (r - gamPre)
      
      ## reset betaPre and gamPre for next lambda if(!is.null(startBeta)) { betaPre=startBeta r=starRes }
      ## if(!is.null(startGam)) { gamPre=startGam }
      
    }  #end iteration for each lambda2 fixed lambda1
    
  }  #end iteration for each lambda1
  
  list(beta = beta, loss = loss, iter = iter, gam = gam)
}

SoftThreshold <- function(z,lambda){
  ifelse(abs(z) <= lambda, 0, ifelse(z >lambda,z-lambda,z+lambda))
}

BIC4PAMLS = function(loss, beta, gam, lambda1, lambda2,criterion = "BIC") {
  # declare and initial
  l1 = length(lambda1)
  l2 = length(lambda2)
  index1 = index2 = 1
  BIC.max= 1e+08
  bicPre = BIC.max
  bicTemp = matrix(0, l1, l2)
  gdf = matrix(0, l1, l2)
  bdf = matrix(0, l1, l2)
  start1 = 1
  end1 = l1
  start2 = 1
  end2 = l2
  n = length(gam[1, 1, ])
  bicTemp2=matrix(0, l1, l2)
  pro=0.5
  for (i in start1:end1) {
    for (j in start2:end2) {
      gdf[i, j] = sum(gam[i, j, ] != 0)
      bdf[i, j] = sum(beta[i, j, ] != 0)
      if (criterion == "AIC") {
        bicTemp[i, j] = log(loss[i, j]/(n)) + (bdf[i, j] + gdf[i, j]) * 2/(n)
      } else {
        # (log(loss/n+a)+log(n)*df/(n))
        bicTemp[i, j] = log(loss[i, j]/(n)) + (bdf[i, j] + gdf[i, j]) * log(n)/(n)
      }
      #compute limit of lambda1 & lambda2
      if(gdf[i,j] >= n * pro || bdf[i,j] + gdf[i,j] >= n){
        bicTemp2[i,j]=-BIC.max
      }else{
        bicTemp2[i,j]=bicTemp[i,j]
        if (bicTemp[i, j] <= bicPre) {
          index1 = i
          index2 = j
          bicPre = bicTemp[i, j]
        }
      }
    }
  }
  bicTemp2 <- ifelse(bicTemp2==-BIC.max,max(bicTemp),bicTemp2)
  res = list(lambda1 = lambda1, lambda2 = lambda2, bdf = bdf, gdf = gdf, bic = bicTemp, bic2=bicTemp2,gam = gam, beta = beta)
  i = index1
  j = index2
  list(lambda1s = lambda1, lambda2s = lambda2, beta = beta[i, j, ], gam = gam[i, j, ], loss = loss[i, j], 
       bdf = bdf[i, j], gdf = gdf[i, j], index1 = i, index2 = j, iter=0,crit1=bicTemp[,j], crit2=bicTemp[i,],res = res)
}

