## This file is for tunning parameter selection

## BIC

BIC_grid = function(wloss, beta, w, lambda1, lambda2, inv = 1, alpha = 0, criterion = "BIC", pro = 0.5) {
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
    bicTemp2=matrix(0, l1, l2)
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
            #compute limit of lambda1 & lambda2
            if(wdf[i,j] >= n * pro || bdf[i,j] + wdf[i,j] >= n){
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
    res = list(lambda1 = lambda1, lambda2 = lambda2, bdf = bdf, wdf = wdf, bic = bicTemp, bic2=bicTemp2,w = w, beta = beta)
    i = index1
    j = index2
    list(lambda1s = lambda1, lambda2s = lambda2, beta = beta[i, j, ], w = w[i, j, ], wloss = wloss[i, j], bdf = bdf[i, 
        j], wdf = wdf[i, j], index1 = i, index2 = j, iter=0,crit1=bicTemp[,j], crit2=bicTemp[i,],res = res)
}


BIC_cross = function(loss, dfw, dfb, n, p, type = "beta", criterion = "BIC", pro = 0.5) {
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
