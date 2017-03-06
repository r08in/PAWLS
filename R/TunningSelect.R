## This file is for tunning parameter selection

## BIC

BICPWLQ = function(wloss, beta, w, lambda1, lambda2, inv = 1, alpha = 0, criterion = "BIC", pro = 0.5) {
    # declare and initial
    l1 = length(lambda1)
    l2 = length(lambda2)
    index1 = index2 = 1
    BIC.max= 1e+08
    bicPre = BIC.max
    bicTemp = matrix(0, l1, l2)
    wdf = matrix(0, l1, l2)
    bdf = matrix(0, l1, l2)
    start1 = 1
    end1 = l1
    start2 = 1
    end2 = l2
    n = length(w[1, 1, ])
    for (i in start1:end1) {
        for (j in start2:end2) {
            wdf[i, j] = sum(w[i, j, ] != 1 + 0)
            bdf[i, j] = sum(beta[i, j, ] != 0 + 0)
            if (criterion == "AIC") {
                bicTemp[i, j] = log(wloss[i, j]/(n)) + (bdf[i, j] + wdf[i, j]) * 2/(n)
            } else {
                # (log(loss/n+a)+log(n)*df/(n))
                bicTemp[i, j] = log(wloss[i, j]/(n)) + (bdf[i, j] + wdf[i, j]) * log(n)/(n)
            }
            if(wdf[i,j] >= n * pro || bdf[i,j] + wdf[i,j] >= n){
              #bicTemp[i,j]=-BIC.max
              next
            }
              
            if (bicTemp[i, j] <= bicPre) {
                index1 = i
                index2 = j
                bicPre = bicTemp[i, j]
            }
        }
    }
    BIC.max=max(bicTemp)
    bicTemp2=matrix(0, l1, l2)
    for (i in start1:end1) {
      for (j in start2:end2) {
        if(wdf[i,j] >= n * pro || bdf[i,j] + wdf[i,j] >= n){
          bicTemp2[i,j]=BIC.max
        }else{
          bicTemp2[i,j]=bicTemp[i,j]
        }
      }
    }
    
    #bicTemp <- ifelse(bicTemp==-BIC.max,max(bicTemp),bicTemp)
    res = list(lambda1 = lambda1, lambda2 = lambda2, bdf = bdf, wdf = wdf, bic = bicTemp, bic2=bicTemp2,w = w, beta = beta)
    i = index1
    j = index2
    list(lambda1s = lambda1, lambda2s = lambda2, beta = beta[i, j, ], w = w[i, j, ], wloss = wloss[i, j], bdf = bdf[i, 
        j], wdf = wdf[i, j], index1 = i, index2 = j, iter=0,crit1=bicTemp[,j], crit2=bicTemp[i,],res = res)
}

BICPWLQ2 = function(wloss, beta, w, lambda1, lambda2, n, inv = 1) {
    # declare and initial
    l1 = length(lambda1)
    l2 = length(lambda2)
    p = dim(beta)[3]
    index1 = index2 = 1
    bicPre = 1e+08
    bicTemp = matrix(0, l1, l2)
    wdf = matrix(0, l1, l2)
    bdf = matrix(0, l1, l2)
    start1 = round(l1 * ((1 - inv)/2)) + 1
    end1 = round(l1 * (0.5 + inv/2))
    start2 = round(l2 * ((1 - inv)/2)) + 1
    end2 = round(l2 * (0.5 + inv/2))
    index1 = l1 - 1
    index2 = BIC(wloss[index1, ], apply(beta[index1, , ] != 0 + 0, 1, sum), n, p, type = "beta")
    index1 = BIC(wloss[, index2], apply(w[, index2, ] != 1 + 0, 1, sum), n, p, type = "w")
    index2 = BIC(wloss[index1, ], apply(beta[index1, , ] != 0 + 0, 1, sum), n, p, type = "beta")
    res = list(lambda1 = lambda1, lambda2 = lambda2, bdf = bdf, wdf = wdf, bic = bicTemp)
    i = index1
    j = index2
    list(lambda1 = lambda1[i], lambda2 = lambda2[j], beta = beta[i, j, ], w = w[i, j, ], wloss = wloss[i, j], bdf = sum(beta[index1, 
        index2, ] != 0 + 0), wdf = sum(w[index1, index2, ] != 1 + 0), index1 = i, index2 = j, res = res)
}

BIC4PAWLS = function(loss, dfw, dfb, n, p, type = "beta", criterion = "BIC", pro = 0.5, blam=NULL) {
    if(n >= p){ # low dimension
      df <- dfb + dfw
      fitLoss <- log(loss/n)
      coeff <- log(n) / n
    }else{ # high dimension
      df <- dfb + dfw
      # fitLoss <- loss/n
      # coeff <- log(n)
      fitLoss <- log(loss/n)
      coeff <- log(n) / n
    }
    if (type == "beta") {
        if (criterion == "AIC") {
            vl = (log(loss/n) + 2 * df/(n))
        } else {
            vl = (fitLoss + coeff * df )
        }
    } else {
        
        if (criterion == "AIC") {
            vl = (log(loss/n) + 2 * df/(n))
        } else {
          vl = (fitLoss + coeff * df)
        }
    }
    crit <- vl
    BIC.max=max(vl)
    for(i in 1: length(vl))
    {
      if(dfw[i] >= n * pro || dfb[i] + dfw[i] >= n)
        crit[i] <- BIC.max
    }
    if(min(crit)==BIC.max){
      crit <- vl
    }
    list(index=which.min(crit),crit=crit)
}

dfs = function(x, beta, w) {
    if (dim(beta)[1] != dim(w)[1]) 
        stop("in df: beta and w should have the same row!")
    require("Matrix")
    m = dim(beta)[1]
    df = rep(0, m)
    b = abs(beta) > 5e-05
    for (i in 1:m) {
        if (sum(b[i, ]) == 0) 
            next
        df[i] = rankMatrix(as.matrix(w[i, ] * x[, b[i, ]]))
    }
    df
}

