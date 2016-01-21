cvreg=function(x, y,penalty1=c("1-w0","log"),penalty2=c("LASSO", "RIDGE", "MCP"),
               lambda1,lambda2,beta0,w0,delta, maxIter,K=5,intercept=TRUE)

{
  n=length(y)
  m=dim(x)[2]
  if(K>n) stop("inside cvreg: k cannot be larger than n")
  size=floor(n/K)
  L1=length(lambda1)
  L2=length(lambda2)
  pe=matrix(0,nrow=L1,ncol=L2)
  rseq=sample(1:n,n)
  
  #interation for each fold
  for(k in 1:K)
  {
    #prepare data
    range=rseq[((k-1)*size+1):(k*size)]
    trainx=x[-range,]
    trainy=y[-range]
    testx=x[range,]
    testy=y[range]
    
    #fit
    res=InnerReg(x=trainx,y=trainy,penalty1=penalty1,penalty2=penalty2,
                 lambda1=lambda1,lambda2=lambda2,
                 beta0=beta0,w0=w0[-range],delta=delta,maxIter=maxIter,intercept=intercept)
    
    #computer pe
    for(l1 in 1:L1)
    {
      for(l2 in 1:L2)
      {
        pe[l1,l2]=pe[l1,l2]+GetRobustPe(X=testx,y=testy,betaHat=res$beta[l1,l2,])
      }
    }

    
  }
  pe=pe/K
  #find best lambda that have testPE smallset
  index1=which.min(pe)%%L1
  index1=ifelse(index1==0,L1,index1)
  index2=ceiling(which.min(pe)/L1)
  res=InnerReg(x=x,y=y,penalty1=penalty1,penalty2=penalty2,lambda1=lambda1[index1],lambda2=lambda2[index2],
               beta0=beta0,w0=w0,delta=delta,maxIter=maxIter,intercept=intercept)
  #return data
  list(beta=as.vector(res$beta),w=res$w,wloss=res$wloss,index1=index1,index2=index2,pe=pe)
}