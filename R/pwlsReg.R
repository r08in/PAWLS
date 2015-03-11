## This functionn is to perform group coordinate descent regression

pwlsreg=function (x,y,lambda=NULL,nlambda=100,w0,delta,maxIter=100,intercept=TRUE,standardize=FALSE,...)
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
  if (nlambda < 2)
    stop("nlambda must be at least 2")
  if (any(is.na(y)) | any(is.na(x))) 
    stop("Missing data (NA's) detected.Take actions to eliminate missing data before passing X and y to gcdreg.")
  w0=ifesle(w0==1,0.99,w0)
  
  #standaredize or not
  if(intercept)
  {
    x=AddIntercept(x)
  }
  std=0
  scale=0
  if(standardize)
  {
    std <- .Call("GroupStandardize", x,y)
    XX <- std[[1]]
    yy <- std[[2]]
    scale <- std[[3]]  
  }
  else
  {
    XX=x
    yy=y
  }
  
  ##setup parameter
  if (missing(lambda)) 
  {
    lambda=SetupParameter1(XX, yy,nlambda,w0,intercept=intercept)
  } 
  else 
  {
    nlambda=length(lambda)
  }
  
  ##Fit  
  res=pwlsReg2(XX, yy,lambda,w0,delta, maxIter,intercept=intercept)
  ##unstandardize 
  if(standardize)
  {
    scale=ifelse(scale==0,0,1/scale)
    res$beta=res$beta*scale
  }
  res
}