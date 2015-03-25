## This functionn is to perform group coordinate descent regression

srcdreg=function (x,y,penalty1="1-w0",penalty2="ADL",lambda1=NULL,lambda2=NULL,nlambda1=100,nlambda2=100,
                 beta0,w0,delta,maxIter=100,intercept=TRUE,standardize=FALSE,updateInitial=TRUE,criterion="BIC",...)#penalty1=c("log","1-w0")
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
  #penalty2 <- match.arg(penalty2)
  #penalty1 <- match.arg(penalty1)
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
  if (missing(lambda1)||missing(lambda2)) 
  {
    lambda=SetupParameter(XX, yy,nlambda1,nlambda2,beta0,w0,intercept=intercept,penalty1=penalty1)
    lambda1=lambda$lambda1
    lambda2=lambda$lambda2
  } 
  else 
  {
    nlambda1=length(lambda1)
    nlambda2=length(lambda2)
  }
  
  ##Fit  
  if(criterion=="BIC")
  {
    res=RCDReg3(XX, yy,penalty1=penalty1,penalty2=penalty2,lambda1,lambda2,beta0,w0,delta, maxIter,intercept=intercept,updateInitial=updateInitial)
  }
  else if(criterion=="CV")
  {
    res=RCDReg4(XX, yy,penalty1=penalty1,penalty2=penalty2,lambda1,lambda2,beta0,w0,delta, maxIter,intercept=intercept,updateInitial=updateInitial)
  }
  ##unstandardize 
  if(standardize)
  {
    scale=ifelse(scale==0,0,1/scale)
    res$beta=res$beta*scale
  }
  res
}