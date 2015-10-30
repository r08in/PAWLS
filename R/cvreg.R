cvreg=function(x, y,penalty1=c("1-w0","log"),penalty2=c("LASSO", "RIDGE", "MCP"),
               lambda1,lambda2,beta0,w0,delta, maxIter,K=10,intercept=TRUE)

{
  n=length(y)
  m=dim(x)[2]
  if(K>n) stop("inside cvreg: K cannot be larger than n")
  size=round(n/K)
  L=length(lambda1)*length(lambda2)
  pe=matrix(nrow=L,ncol=K,0) #robust
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
    res=InnerReg(trainx,trainy,penalty1=penalty1,penalty2=penalty2,lambda1,lambda2,
             beta0=rep(1,m),w0=rep(1,n)[-range],delta,maxIter,intercept=intercept)
    beta=matrix(res$beta,L,m)
    
    #computer pe
    for(i in 1:L)
    {
      pe[i,k]=GetRobustPe(X=testx,y=testy,betaHat=beta[i,])
    }
    
  }
  #find best lambda that have testPE smallset
  pe=apply(pe,1,sum)/K
  index=which.min(pe)
  if(length(lambda1)==1)
  {
    lambda2=lambda2[index]
  }
  else 
  {
    lambda1=lambda1[index]
  }
  res=InnerReg(x,y,penalty1=penalty1,penalty2=penalty2,lambda1,lambda2,
           beta0,w0,delta,maxIter,intercept=intercept)
  #return data
  list(beta=res$beta,w=res$w,wloss=res$wloss,index=index,pe=pe)
}