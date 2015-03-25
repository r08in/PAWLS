RCDReg4=function(x, y,penalty1=c("log","1-w0"),penalty2=c("MCP", "SCAD", "ADL"),lambda1,lambda2,beta0,w0,delta, 
                 maxIter,intercept=TRUE,updateInitial=TRUE,K=5)
{
  
  L1=length(lambda1)
  L2=length(lambda2)
  m=dim(x)[2]
  n=length(y)
  pre1=pre2=0
  lmaxIter=10
  
  #initial lambda
  index2=L2/2
  res=cvreg(x,y,penalty1=penalty1,penalty2=penalty2,lambda1,lambda2[index2],
            beta0,w0,delta,maxIter,intercept=intercept,K=K) #fix lambda2 for lambda1
  index1=res$index
  res=cvreg(x,y,penalty1=penalty1,penalty2=penalty2,lambda1[index1],lambda2,
           beta0,w0,delta,maxIter,intercept=intercept,K=K) #fix lambda1 for lambda2
  index2=res$index
  ##loop to estimate and find the best
  iter=0
  while((pre1!=index1||pre2!=index2)&&(iter<lmaxIter))
  {
    iter=iter+1
    pre1=index1
    pre2=index2
    
    if(updateInitial)
    {
      beta0=res$beta    
      w0=res$w
      beta0=SetBeta0(beta0)
      w0=ifelse(w0==1,0.99,w0)
    }    
    res=cvreg(x,y,penalty1=penalty1,penalty2=penalty2,lambda1,lambda2[index2],
              beta0,w0,delta,maxIter,intercept=intercept,K=K) #fix lambda2 for lambda 1
    index1=res$index
    if(updateInitial)
    {
      beta0=res$beta    
      w0=res$w
      beta0=SetBeta0(beta0)
      w0=ifelse(w0==1,0.99,w0)
    }    
    res=cvreg(x,y,penalty1=penalty1,penalty2=penalty2,lambda1[index1],lambda2,
                 beta0,w0,delta,maxIter,intercept=intercept) #fix lambda1 for lambda2
    index2=res$index
    if(pre2==index2&pre1==index1) break;
  }
  
  if(updateInitial)
  {
    beta0=res$beta    
    w0=res$w
    beta0=SetBeta0(beta0)
    w0=ifelse(w0==1,0.99,w0)
  }  
  #test##
  res=InnerReg(x,y,penalty1=penalty1,penalty2=penalty2,lambda1[index1],lambda2[index2],beta0,w0,delta,maxIter,intercept=intercept) #fix lambda2 and lambda2
  #return 
  i=index1
  j=index2
  list(lambda1=lambda1[i],lambda2=lambda2[j],
       beta=as.vector(res$beta),w=as.vector(res$w),
       wloss=res$wloss,bdf=sum(res$beta!=0+0),wdf=sum(res$w!=1+0),
       index1=i,index2=j,iter=iter)
  
}