## This file is for tunning parameter selection

## BIC
BIC_grid = function(wloss, beta, w) {
  
  l1 <- dim(beta)[1]
  l2 <-  dim(beta)[2]
  index1 <- 1
  index2 <- 1
  BIC.max <- 1e+08
  bicPre = BIC.max
  bicTemp <- matrix(0, l1, l2)
  bicTemp2<- matrix(0, l1, l2)
  wdf <- matrix(0, l1, l2)
  bdf <- matrix(0, l1, l2)
  n <- length(w[1, 1, ])
  pro <- 0.5
  
  for (i in 1 : l1) {
    for (j in 1 : l2) {
      wdf[i, j] <- sum(w[i, j, ] != 1 + 0)
      bdf[i, j] <- sum(beta[i, j, ] != 0 + 0)
      bicTemp[i, j] <- log(wloss[i, j]/(n)) + (bdf[i, j] + wdf[i, j]) * log(n)/(n)
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
  res = list(beta = beta[index1, index2, ], w = w[index1,index2,], raw.bic = bicTemp, 
             bic=bicTemp2,index1=index1,index2=index2)
}

BIC_cross = function(loss, dfw, dfb, n ) {
  pro <- 0.5
  df <- dfb + dfw
  fitLoss <- log(loss/n)
  coeff <- log(n) / n
  vl = (fitLoss + coeff * df )
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
