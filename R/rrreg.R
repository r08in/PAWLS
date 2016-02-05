## This functionn is to perform group coordinate descent regression
#penalty1=c("log","1-w0")
#initial=
rrreg=function (x,y,penalty1=c("1-w0","log"),penalty2=c("RIDGE"),
                  lambda1=NULL,lambda2=NULL,nlambda1=30,nlambda2=30,
                  beta0=rep(1,p),w0=rep(0,n),initial=c("RRMM"),
                  delta=0.000001,maxIter=1000,matlab=NULL,
                  intercept=TRUE,standardize=FALSE,
                  criterion=c("CV","BIC","AIC"),dout=FALSE,...)
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
  penalty2 <- match.arg(penalty2)
  initial <- match.arg(initial)
  criterion <- match.arg(criterion)
  
  if (nlambda1 < 1||nlambda2<1) 
    stop("nlambda must be at least 1")
  if(!is.null(lambda1)) 
    nlambda1=length(lambda1)
  if(!is.null(lambda2)) 
    nlambda2=length(lambda2)
  
  if (any(is.na(y)) | any(is.na(x))) 
    stop("Missing data (NA's) detected.Take actions to eliminate missing data before passing X and y to gcdreg.")
  

  #intercept
#   if(intercept)
#   {
#     x=AddIntercept(x)
#   }
  
  #initial
  n=length(y)
  p=dim(x)[2]
  
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
#   if (is.null(lambda1))
#   {
#     if(initial=="RRMM")
#     {
#       betaInit=GetBetaInit(XX,yy)
#     }
#     lambda1=InitLambda1(XX,yy,betaInit=betaInit)
#   } 
#   if (is.null(lambda2))
#   {
#     lambda2=seq(1,0.001,length.out=100)
#   } 
  ##Fit  
  
  if(criterion=="BIC"||criterion=="AIC")
  {
  
    
  }
  else if(criterion=="CV")
  {
    #res=cvreg(XX, yy,penalty1=penalty1,penalty2=penalty2,lambda1,lambda2,beta0,w0,delta, maxIter,intercept=intercept)
    res=cvreg2(x=XX, y=yy,lambda1=lambda1,lambda2=lambda2,beta0=beta0,w0=w0,matlab=matlab,dout=dout)
    #plot
#     nlambda1=length(lambda1)
#     nlambda2=length(lambda2)
#     pe=matrix(0,nrow=(nlambda1*nlambda2),ncol=3)
#     for(l1 in 1:nlambda1)
#     {
#       for(l2 in 1:nlambda2)
#       {
#         pe[(l2-1)*nlambda1+l1,1]=lambda1[l1]
#         pe[(l2-1)*nlambda1+l1,2]=lambda2[l2]
#         pe[(l2-1)*nlambda1+l1,3]=res$pe[l1,l2]
#       }
#     }
    #scatter3d(pe[,1],pe[,3],pe[,2])
    #\plot3d(pe[,1],pe[,2],pe[,3])
  }
  ##unstandardize 
  
  if(standardize)
  {
    scale=ifelse(scale==0,0,1/scale)
    res$beta=res$beta*scale
  }
  res
}