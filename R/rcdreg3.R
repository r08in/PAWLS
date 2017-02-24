# without search the whole grid of lambda
RCDReg3 = function(x, y, penalty1 = c("1-w0", "log"), penalty2 = c("LASSO", "RIDGE", "MCP"), lambda1, lambda2, beta0, 
    w0, delta, maxIter, intercept = TRUE, updateInitialTimes = 0, criterion = "BIC", startBeta = NULL, startW = NULL) {
    penalty1 <- match.arg(penalty1)
    penalty2 <- match.arg(penalty2)
    
    L1 = length(lambda1)
    L2 = length(lambda2)
    m = dim(x)[2]
    n = length(y)
    pre1 = pre2 = 0
    lmaxIter = 10
    beta0 = SetBeta0(beta0)
    w0 = ifelse(w0 == 1, 0.99, w0)
    
    # initial lambda
    index2 = round(L2/2)
    res = InnerReg(x, y, penalty1 = penalty1, penalty2 = penalty2, lambda1, lambda2[index2], beta0, w0, delta, maxIter, 
        intercept = intercept, startBeta = startBeta, startW = startW)  #fix lambda2
    bic1 = BIC4PAWLS(as.vector(res$wloss), apply(matrix(res$w, L1, n) != 1 + 0, 1, sum), 
                     apply(matrix(res$beta, L1, m) != 0 + 0, 1, sum), n, m, type = "w", criterion = criterion, 
                     x = x, ws = matrix(res$w, L1, n), bs = matrix(res$beta, L1, m)) 
    index1 =bic1$index  #find best lambda1
    crit1 = bic1$crit

    res = InnerReg(x, y, penalty1 = penalty1, penalty2 = penalty2, lambda1[index1], lambda2, beta0, w0, delta, maxIter, 
        intercept = intercept, startBeta = startBeta, startW = startW)  #fix lambda1
    bic2 = BIC4PAWLS(as.vector(res$wloss), apply(matrix(res$w, L2, n) != 1 + 0, 1, sum),
                     apply(matrix(res$beta,L2, m) != 0 + 0, 1, sum), n, m, 
                     criterion = criterion, x = x, ws = matrix(res$w, L2, n), bs = matrix(res$beta, L2, m))
    index2 = bic2$index  #find best lambda2
    crit2 = bic2$crit
    ## loop to estimate and find the best
    iter = 0
    
    while ((pre1 != index1 || pre2 != index2) && (iter < lmaxIter)) {
        iter = iter + 1
        pre1 = index1
        pre2 = index2
        
        if (iter <= updateInitialTimes) {
            beta0 = matrix(res$beta, L2, m)[index2, ]
            w0 = matrix(res$w, L2, n)[index2, ]
            beta0 = SetBeta0(beta0)
            w0 = ifelse(w0 == 1, 0.99, w0)
        }
        resw = InnerReg(x, y, penalty1 = penalty1, penalty2 = penalty2, lambda1, lambda2[index2], beta0, w0, delta, 
            maxIter, intercept = intercept, startBeta = startBeta, startW = startW)  #fix lambda2
        bic1 = BIC4PAWLS(as.vector(resw$wloss), apply(matrix(resw$w, L1, n) != 1 + 0, 1, sum), 
                         apply(matrix(resw$beta, L1, m) != 0 + 0, 1, sum), n, m, type = "w", criterion = criterion, 
                         x = x, ws = matrix(resw$w, L1, n), bs = matrix(resw$beta, L1, m)) 
        index1 =bic1$index  #find best lambda1
        crit1 = bic1$crit
        
        if (iter <= updateInitialTimes) {
            beta0 = matrix(resw$beta, L1, m)[index1, ]
            w0 = matrix(resw$w, L1, n)[index1, ]
            beta0 = SetBeta0(beta0)
            w0 = ifelse(w0 == 1, 0.99, w0)
        }
        res = InnerReg(x, y, penalty1 = penalty1, penalty2 = penalty2, lambda1[index1], lambda2, beta0, w0, delta, 
            maxIter, intercept = intercept, startBeta = startBeta, startW = startW)  #fix lambda1
        bic2 = BIC4PAWLS(as.vector(res$wloss), apply(matrix(res$w, L2, n) != 1 + 0, 1, sum),
                         apply(matrix(res$beta,L2, m) != 0 + 0, 1, sum), n, m, 
                         criterion = criterion, x = x, ws = matrix(res$w, L2, n), bs = matrix(res$beta, L2, m))
        index2 = bic2$index  #find best lambda2
        crit2 = bic2$crit
        if (pre2 == index2 && pre1 == index1) 
        {
          break
        }
            
    }
    
    list(lambda1 = lambda1[index1], lambda2 = lambda2[index2], beta = as.vector(res$beta[1, index2, ]), w = as.vector(res$w[1, 
        index2, ]), wloss = res$wloss[1, index2], bdf = sum(res$beta[1, index2, ] != 0 + 0), wdf = sum(res$w[1, index2, 
        ] != 1 + 0), index1 = index1, index2 = index2, iter = iter, w0 = w0, beta0 = beta0, lambda1s = lambda1, lambda2s = lambda2, 
        betas = matrix(res$beta, L2, m), ws = matrix(resw$w, L1, n), crit1=crit1, crit2=crit2)
    
}

InnerReg = function(x, y, penalty1 = "1-w0", penalty2 = c("LASSO", "RIDGE", "MCP"), lambda1, lambda2, beta0, w0, 
    delta, maxIter, intercept = TRUE, startBeta = NULL, startW = NULL) {
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
        startW = rep(1, n)
    }
    res <- .Call("INNERREG", x, y, penalty1, penalty2, lambda1, lambda2, beta0, w0, delta, maxIter, ifelse(intercept, 
        1, 0), startBeta = startBeta, startW = startW)
    
    res = list(beta = array(res[[1]], dim = c(L1, L2, m)), w = array(res[[2]], dim = c(L1, L2, n)), wloss = array(res[[3]], 
        dim = c(L1, L2)), loss = array(res[[4]], dim = c(L1, L2)), iter = array(res[[5]], dim = c(L1, L2)))
    
    res
    
}
InnerReg2 = function(x, y, penalty1 = "1-w0", penalty2 = "ADL", lambda1, lambda2, beta0, w0, delta, maxIter, intercept = TRUE, 
    startBeta = NULL, startW = NULL) {
    
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
    w = array(1, dim = c(L1, L2, n))
    loss = matrix(ncol = L2, nrow = L1, 0)
    wloss = matrix(ncol = L2, nrow = L1, 0)
    iter = matrix(ncol = L2, nrow = L1, 0)
    
    betaPre = rep(0, m)
    wPre = rep(1, n)
    if (!is.null(startBeta)) {
        betaPre = startBeta
    }
    if (!is.null(startW)) {
        wPre = startW
    }
    starRes = y - x %*% betaPre
    r = starRes
    ## iteration for each lamda1
    for (l1 in lstart1:L1) {
        ## initial
        
        shift = rep(0, m + n)
        c = rep(0, m)
        loss[l1, 1] = t(y) %*% y  ##initial loss[l1,1]
        if (penalty1 == "log") {
            lam1 = sqrt((lambda1[l1]/abs(log(w0))) * n)  #init sqrt(lambda1/abs(log(w0))n)
        } else {
            lam1 = (lambda1[l1]/abs(1 - w0)) * n
        }
        
        ## iteration for each lamda2
        for (l2 in lstart2:L2) {
            lam2 = lambda2[l2]/abs(beta0)
            if (intercept) {
                lam2[1] = 0
            }
            ## iteration for all covariates
            while (iter[l1, l2] < maxIter) {
                iter[l1, l2] = iter[l1, l2] + 1
                
                # calculate coefficient c
                c = apply((x * wPre)^2, 2, sum)/n
                
                ## iteration for each beta
                for (j in 1:m) {
                  ## (1)calculate zj
                  zj = t(x[, j] * wPre^2) %*% r/n + c[j] * betaPre[j]
                  ## (2)update betaj
                  if (penalty2 == "LASSO") {
                    beta[l1, l2, j] = UpdateBeta(zj, lam2[j], c[j])
                    # beta[[l1]][l2,j]=UpdateBeta(zj,lam2[j]*sqrt(c[j]),c[j]) beta[[l1]][l2,j]=zj
                  }
                  
                  ## (3)update r
                  shift[j] = beta[l1, l2, j] - betaPre[j]
                  r = r - x[, j] * shift[j]
                }
                
                ## update w
                if (penalty1 == "log") {
                  absr = abs(r)
                  w[l1, l2, ] = ifelse(absr > lam1, lam1/absr, 1)
                } else {
                  # 1-w0
                  sqr = r^2
                  w[l1, l2, ] = ifelse(sqr > lam1, lam1/sqr, 1)
                }
                shift[(m + 1):(m + n)] = w[l1, l2, ] - wPre
                
                ## update betaPre and wPre for next iteration
                betaPre = beta[l1, l2, ]
                wPre = w[l1, l2, ]
                ## Check for convergence
                if (t(shift) %*% shift < delta) {
                  break
                }
                
            }  #end for the inner loop
            
            ## compute square of loss
            loss[l1, l2] = t(r) %*% r
            wloss[l1, l2] = t(r * wPre) %*% (r * wPre)
            
            ## reset betaPre and wPre for next lambda if(!is.null(startBeta)) { betaPre=startBeta r=starRes }
            ## if(!is.null(startW)) { wPre=startW }
            
        }  #end iteration for each lambda2 fixed lambda1
        
    }  #end iteration for each lambda1
    
    list(beta = beta, loss = loss, iter = iter, w = w, wloss = wloss)
}


