## This functionn is to perform robust group coordinate descent regression

rgcdreg=function (x,y,groupInfo,penalty=c("MCP", "SCAD", "ADL"),gamma,lambda1=NULL,lambda2=NULL,nlambda1=100,nlambda2=100,
                  w0,delta,maxIter=100,...)
{
  ##error checking
  if (class(x) != "matrix") 
  {
    tmp <- try(x <- as.matrix(x), silent=TRUE)
    if (class(tmp)[1] == "try-error") 
      stop("x must be a matrix or able to be coerced to a matrix")
  }
  if (class(y) != "numeric") 
  {
    tmp <- try(y <- as.numeric(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") 
      stop("y must numeric or able to be coerced to numeric")
  }
  penalty <- match.arg(penalty)
  if (gamma <= 1 & penalty=="MCP") 
    stop("gamma must be greater than 1 for the MC penalty.")
  if (gamma <= 2 & penalty=="SCAD") 
    stop("gamma must be greater than 2 for the SCAD penalty.")
  if (nlambda1 < 2||nlambda2<2) 
    stop("nlambda must be at least 2")
  if (any(is.na(y)) | any(is.na(x))) 
    stop("Missing data (NA's) detected.Take actions to eliminate missing data before passing X and y to gcdreg.")
  
  ##group standardize
  std <- .Call("GroupStandardize", x,y)
  XX <- std[[1]]
  yy <- std[[2]]
  scale <- std[[3]]  
  
  
  ##setup parameter
  if (missing(lambda1)||missing(lambda2)) 
  {
    lambda=SetupGroupParameter(XX,yy,nlambda1,nlambda2,w0,groupInfo)
    if(missing(lambda1)) lambda1=lambda$lambda1
    if(missing(lambda2)) lambda2=lambda$lambda2
    
  } 
  nlambda1=length(lambda1)
  nlambda2=length(lambda2)
  
  ##Fit
  res=RGCDReg2(XX,yy,groupInfo,penalty,gamma,lambda1,lambda2,w0,delta, maxIter)
  m<-ncol(XX)
  n<-nrow(XX)
  
  
  ##unstandardize
  scale=ifelse(scale==0,0,1/scale)
  for(i in 1:length(res$beta))
  {
    res$beta[[i]]=t(t(res$beta[[i]])*(scale))
  }
  
  ##output
  val <- structure(list(beta = res$beta,
                        w=res$w,
                        iter = res$iter,
                        lambda1 = lambda1,
                        lambda2=lambda2,
                        penalty = penalty,
                        loss = res$loss,
                        wloss= res$wloss,
                        n = n),
                   class = "rgcdreg")
  val
}