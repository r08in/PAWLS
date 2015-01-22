SetupParameter = function(x,y,nlambda1,nlambda2,beta0,w0) 
{
  #set lambda1
  lambda1Max=max(y^2*abs(log(w0))/n)
  lambda1=logSeq(lambda1Max,lambda1Max*0.0001,nlambda1)
  #set lambda2
  lambda2Max=max(abs(t(x)%*%y/n)*abs(beta0)) # max |betaj|*|xj'y/n|
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
