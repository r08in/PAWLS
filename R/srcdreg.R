## This functionn is to perform group coordinate descent regression

srcdreg=function (x,y,penalty=c("MCP", "SCAD", "ADL"),lambda1=NULL,lambda2=NULL,nlambda1=100,nlambda2=100,
                 beta0,w0,delta,maxIter=100,...)
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
  if (nlambda1 < 2||nlambda2<2) 
    stop("nlambda must be at least 2")
  if (any(is.na(y)) | any(is.na(x))) 
    stop("Missing data (NA's) detected.Take actions to eliminate missing data before passing X and y to gcdreg.")
  
  ##standardize
  #std <- .Call("standardize", x,y)
  #XX <- std[[1]]
  #center <- std[[2]]
  #scale <- std[[3]]
  #yy <- y - mean(y)
  
  std <- .Call("GroupStandardize", x,y)
  XX <- std[[1]]
  yy <- std[[2]]
  scale <- std[[3]]  
  
  
  ##setup parameter
  if (missing(lambda1)||missing(lambda2)) 
  {
    lambda=SetupParameter(XX, yy,nlambda1,nlambda2,beta0,w0)
    lambda1=lambda$lambda1
    lambda2=lambda$lambda2
  } 
  else 
  {
    nlambda1=length(lambda1)
    nlambda2=length(lambda2)
  }
  
  ##Fit  
  res=RCDReg3(XX, yy,penalty,lambda1,lambda2,beta0,w0,delta, maxIter)
  ##unstandardize 
  scale=ifelse(scale==0,0,1/scale)
  res$beta=res$beta*scale
  res
}