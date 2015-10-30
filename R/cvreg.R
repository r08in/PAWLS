cvreg=function(x, y,penalty1="1-w0",penalty2="ADL",lambda1,lambda2,beta0,w0,delta, maxIter,
               standardiz=FALSE,intercept=TRUE,K=10)
{
  n=length(y)
  m=dim(x)[2]
  if(K>n) stop("inside cvreg: K cannot be larger than n")
  size=round(n/K)
  L=length(lambda1)*length(lambda2)
  pe=matrix(nrow=L,ncol=K,0)
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
             beta0,w0[-range],delta,maxIter,intercept=intercept)
    beta=matrix(res$beta,L,m)
    
    #computer pe
    for(i in 1:L)
    {
      pe[i,k]=sum((testy-testx%*%beta[i,])^2)
    }
    
  }
  #find best lambda that have testPE smallset
  index=which.min(apply(pe,1,sum))
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
  list(index=index,beta=res$beta,w=res$w,wloss=res$wloss)
}