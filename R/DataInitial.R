##initial

InitBeta=function(x,y,method="LAD",nlambda=300,selector="BIC")
{
  res=NULL
  require("flare")
  if(method=="LAD")
  {    
    res=slim(x,y,method="lq",q=1,nlambda=nlambda,verbose=FALSE)
  }
  else if(method=="LASSO")
  {
    res=slim(x,y,method="lasso",nlambda=nlambda,verbose=FALSE)
  }
  
  #get loss
  r=y-x%*%res$beta
  loss=apply(r*r,2,sum)
  
  if(selector=="BIC")
  {
    out=BICSelect2(loss,length(y),t(res$beta),res$lambda)
  }
  return(out$beta)
}

UpdateWeight=function(x,y,beta,lambda1=NULL)
{
  n=length(y)
  r=abs(y-x%*%beta)
  if(is.null(lambda1))
  {
    z=sqrt(sum(r*r)/n)
  }
  else
  {
    z=sqrt(lambda1*n)
  }
  w=ifelse(r>z,z/r,1)
}

InitParam=function(x,y,method="LAD",nlambda=300,selector="BIC",lambda1=NULL)
{
  beta=InitBeta(x,y,method=method,nlambda=nlambda,selector=selector)
  weight=UpdateWeight(x,y,beta=beta,lambda1=lambda1)
  return(list(beta=beta,weight=weight))
}