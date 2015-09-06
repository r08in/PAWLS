## This functionn is to perform group coordinate descent regression 
## use for adaptive lasso now

rcdreg=function (x,y,penalty1="1-w0",penalty2="ADL",lambda1=NULL,lambda2=NULL,
                 nlambda1=100,nlambda2=100,
                  beta0,w0,delta,maxIter=100,
                 intercept=TRUE,standardize=FALSE,
                 updateInitial=TRUE,criterion="BIC",...)#penalty1=c("log","1-w0")
{
  ###error checking
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
  res=RCDReg2(XX, yy,penalty1=penalty1,penalty2=penalty2,lambda1,lambda2,beta0,w0,delta, maxIter,intercept=intercept)
  m<-ncol(XX)
  n<-nrow(XX)

  
  ##unstandardize
  ##unstandardize 
  if(standardize)
  {
    scale=ifelse(scale==0,0,1/scale)
    for(i in 1:dim(res$beta)[1])
    {
      res$beta[i,,]=t(t(res$beta[i,,])*(scale))
    }
  }
 
  
  ##output
  BICPWLQ(res$wloss,res$beta,res$w,lambda1,lambda2,n,inv=1)
  #val <- structure(list(beta = res$beta,
                       # w=res$w,
                        #iter = res$iter,
                        #lambda1 = lambda1,
                        #lambda2=lambda2,
                        #loss = res$loss,
                        #wloss= res$wloss,
                       # n = n),
                   #class = "rcdreg")
  #val
}