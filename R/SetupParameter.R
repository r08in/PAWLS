SetupParameter = function(x,y,nlambda1,nlambda2,beta0,w0,intercept=TRUE,alpha=0.1,penalty1="1-w0") 
{
  #set lambda2
  if(intercept)
  {
    #set lambda1
    if(penalty1=="log")
    {
      l1=(y-mean(y))^2*abs(log(w0))/n
      
    }
    else #penalty1="1-w0"
    {
      l1=(y-mean(y))^2*abs(1-w0)/n
    }
    #
    num=round(length(y)*alpha)
    lambda1Max=l1[order(l1,decreasing=TRUE)[num+1]]
    
    lambda2Max=max(abs(t(x[,-1])%*%(y-mean(y))/n)*abs(beta0[-1])) # max |betaj|*|xj'y/n|
  }
  else
  {
    if(penalty1=="log")
    {
      l1=y^2*abs(log(w0))/n
    }
    else #penalty1="1-w0"
    {
      l1=y^2*abs(1-w0)/n
    }
    num=round(length(y)*alpha)
    lambda1Max=l1[order(l1,decreasing=TRUE)[1]]
    lambda2Max=max(abs(t(x)%*%y/n)*abs(beta0)) # max |betaj|*|xj'y/n|
  }
 
  lambda1=logSeq(lambda1Max,lambda1Max*0.01,nlambda1)
  lambda2=logSeq(lambda2Max,0,nlambda2)
  return(list(lambda1=lambda1,lambda2=lambda2))
}

logSeq=function(smax,smin,n)
{
  smin+(smax-smin)/(log(n))*(log(n)-log(1:n))
}


SetupGroupParameter = function(x,y,nlambda1,nlambda2,w0,groupInfo) 
{
  #set lambda1
  lambda1Max=max(y^2*abs(log(w0))/n)
  lambda1=logSeq(lambda1Max,0,nlambda1)
  
  #set lambda2
  lambda2Max=(.Call("MaxProduct",x, y, groupInfo)/sqrt(min(groupInfo)))  
  lambda2=logSeq(lambda2Max,0,nlambda2)
  return(list(lambda1=lambda1,lambda2=lambda2))
}
