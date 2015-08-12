## This functionn is to perform group coordinate descent regression
#penalty1=c("log","1-w0")
#initial=
srcdreg=function (x,y,penalty1=c("1-w0","log"),penalty2="ADL",
                  lambda1=NULL,lambda2=NULL,nlambda1=50,nlambda2=100,
                  beta0=NULL,w0=NULL,initial=c("norm","LTS"),
                  delta=0.000001,maxIter=1000,
                  intercept=TRUE,standardize=FALSE,
                  updateInitialTimes=0,criterion=c("BIC","CV"),...)
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
  
  penalty1 <- match.arg(penalty1)
  initial <- match.arg(initial)
  criterion <- match.arg(criterion)
  
  if(!is.null(lambda1)) 
    nlambda1=length(lambda1)
  if(!is.null(lambda2)) 
    nlambda2=length(lambda2)
  if (nlambda1 < 2||nlambda2<2) 
    stop("nlambda must be at least 2")
  
  if (any(is.na(y)) | any(is.na(x))) 
    stop("Missing data (NA's) detected.Take actions to eliminate missing data before passing X and y to gcdreg.")
  
  #initial
  if(initial=="LTS")
  {
    require(robustHD)
    init=sparseLTS(x,y,intercept=intercept)
    beta0=SetBeta0(init$coefficients)
    w0=ifelse(init$wt==1,0.99,0.01)
    #w0=UpdateWeight(x,y,beta0)
    #w0=ifelse(as.vector(w0)==1,0.99,w0)
  }
  else if(initial=="norm")
  {
    if(is.null(beta0))
    {
      beta0=rep(1,p)
      if(intercept) beta0=c(1,beta0)
    }
    if(is.null(w0))
      w0=rep(0.99,n)
  }
  
  #intercept
  if(intercept)
  {
    x=AddIntercept(x)
  }
  
  #sandardize
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
    if(is.null(lambda1))
      lambda1=lambda$lambda1
    if(is.null(lambda2))
      lambda2=lambda$lambda2
  } 
  
  ##Fit  
 
  if(criterion=="BIC")
  {
    res=RCDReg3(XX, yy,penalty1=penalty1,penalty2=penalty2,lambda1,lambda2,beta0,w0,delta, maxIter,intercept=intercept,updateInitialTimes=updateInitialTimes)
  }
  else if(criterion=="CV")
  {
    res=RCDReg4(XX, yy,penalty1=penalty1,penalty2=penalty2,lambda1,lambda2,beta0,w0,delta, maxIter,intercept=intercept)
  }
  ##unstandardize 
  if(standardize)
  {
    scale=ifelse(scale==0,0,1/scale)
    res$beta=res$beta*scale
  }
  res
}