InitLambda1=function(X,y,betaInit,coeff=1,pro=0.1)
{
  num=round(length(y)*(1-pro))
  p=dim(X)[2]
  if(p>=num)
    stop("number of valid observations can't be less than p!")
  X=AddIntercept(X)
  resid=y-X%*%betaInit
  resid2=resid[order(resid,decreasing=FALSE)[num]]
  sigmaHat=sum(resid2^2)/(num-p)
  coeff*sigmaHat
}

GetRobustPe=function(X,y,betaHat,pro=0.1)
{
  num=round(length(y)*(1-pro))
  p=dim(X)[2]
  resid=y-X%*%betaHat
  resid2=resid[order(resid,decreasing=FALSE)[num]]
  sum(resid2^2)/num
}

GetBetaInit=function(X,y,type=c("RRMM"))
{
  if(type=="RRMM")
  {
    setVariable(matlab, X=X)
    setVariable(matlab, y=y)
    evaluate(matlab,"[beta resid edf lamin]=RobRidge(X,y)")
    betaRRMM=getVariable(matlab, "beta")
    as.vector(betaRRMM$beta)
  }
  
}
