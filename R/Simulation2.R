simulate2=function(L,n,beta=NULL,model=c("A","B","C","D"),p=NULL,method="PAWLS",matlab=NULL,seed=2014,useDataFile=FALSE,
                  standardize=FALSE,penalty1="1-w0",penalty2="RIDGE",updateInitial=FALSE,criterion="BIC",intercept=TRUE,initial="uniform",
                  range="cross",type=c("Lasso","Ridge"),lambda1=NULL,lambda2=NULL)
{
  ptm=proc.time()
  if(!is.null(seed))
  {
    set.seed(seed)
  }
  p=length(beta)
  if(useDataFile)
  {
    f=paste("data\\",model,n,"X",p,".rda",sep="")
    load(f)
    beta=data[[1]]$beta
    n=length(data[[1]]$y)
  }
  pe=rep(0,L)
  w=matrix(0,nrow=L,ncol=n)
  b=matrix(0,nrow=L,ncol=p+1)
  #out=GenerateDataByModel(n=n,beta=beta,model=model)
  #evaluate(matlab,"rng(2015);")
  for(i in 1:L)
  {
    if(i==16)
    {
      test=1+1
    }
    if(useDataFile)
    {
      out=data[[i]]
    }
    else
    {
      out=GenerateDataByModel(n=n,beta=beta,model=model,dataType=type)
    }
    
    ####### for Ridge#######
    if(method=="RRMM")
    {
      setVariable(matlab, X=out$x)
      setVariable(matlab, y=out$y)
      evaluate(matlab,"[beta resid edf lamin]=RobRidge(X,y)")
      betaRRMM=getVariable(matlab, "beta")
      b[i,]=c(betaRRMM$beta[p+1],betaRRMM$beta[1:p])
    }
    else if(method=="RRREG")
    {
      res=rrreg(x=out$x,y=out$y,lambda1=lambda1,lambda2=lambda2,intercept=intercept)
      b[i,]=res$beta
      w[i,]=as.vector(res$w)
    }
    
    #pe
    x=AddIntercept(out$x)
    pe[i]=GetRobustPe(x=x,y=out$y,betaHat=b[i,] )
  }

  
  time=(proc.time()-ptm)[1]
  #return
  list(pe=pe,time=time,w=w,beta=b)
}