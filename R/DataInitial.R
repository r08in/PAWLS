##initial

BetaInit=function(x,y,method="LAD",nlambda=300,selector="BIC")
{
  res=NULL
  require("flare")
  if(method=="LAD")
  {    
    res=slim(x,y,method="lq",q=1,nlambda=nlambda)
  }
  else if(method=="LASSO")
  {
    res=slim(x,y,method="lasso",nlambda=nlambda)
  }
  
  #get loss
  r=y-x%*%res$beta
  loss=apply(r*r,2,sum)
  
  if(selector=="BIC")
  {
    out=BICSelect2(loss,length(y),t(res$beta),res$lambda)
  }
  return(out)
}