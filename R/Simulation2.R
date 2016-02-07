simulate2=function(L,n,beta=NULL,model=c("A","B","C","D"),p=NULL,method="PAWLS",matlab=NULL,seed=2014,useDataFile=FALSE,
                  standardize=FALSE,penalty1="1-w0",penalty2="RIDGE",updateInitial=FALSE,criterion="BIC",intercept=TRUE,initial="uniform",
                  range="cross",type=c("Lasso","Ridge"),lambda1=NULL,lambda2=NULL,r=0.9,errorSigma=1,pro=0.1,dout=FALSE)
{
  ptm=proc.time()
  if(!is.null(seed))
  {
    set.seed(seed)
  }
  p=length(beta)
  
  #prepare matlab
  if(initial=="RRMM"||method=="RRMM")
  {
#     matLabDir="D:\\matlab\\RRMM"
#     matlab=PrepareMatlab(matLabDir)
  }
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
  index1=rep(0,L)
  index2=rep(0,L)
  iter=rep(0,L)
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
      out=GenerateDataByModel(n=n,beta=beta,model=model,r=r,errorSigma=errorSigma,dataType=type,pro=pro)
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
      if(initial=="RRMM")
      {
        setVariable(matlab, X=out$x)
        setVariable(matlab, y=out$y)
        evaluate(matlab,"[beta resid edf lamin]=RobRidge(X,y)") #with intercept
        betaRRMM=getVariable(matlab, "beta")
        beta0=c(betaRRMM$beta[p+1],betaRRMM$beta[1:p])
        #beta0=c(-0.0416693,2.5240519,1.7076847,0.9545043,1.1682181,0.9077817,1.3099283,0.6452162,0.2643470)
        sqr=(out$y-AddIntercept(out$x)%*%beta0)^2
        lam1=1 #lambda1*n
        w0=ifelse(sqr>lam1,lam1/sqr,0.99)
      }
      res=rrreg(x=out$x,y=out$y,lambda1=lambda1,lambda2=lambda2,beta0=beta0,w0=w0,intercept=intercept,matlab=matlab,dout=dout)
      b[i,]=res$beta
      w[i,]=as.vector(res$w)
      index1[i]=res$index1
      index2[i]=res$index2
      iter[i]=res$iter
    }
    else if(method=="ARRREG")
    {
      init=rrreg(x=out$x,y=out$y,lambda1=seq(from=1,to=2,length.out=100),
                lambda2=seq(from=0.005,to=0.5,length.out=100),intercept=intercept)
      b[i,]=res$beta
      w[i,]=as.vector(res$w)
      index1[i]=res$index1
      index2[i]=res$index2
      
    }
  }

  #MSE for beta
  #SE
  se=apply((t(b)-c(0,beta))^2,2,sum)
  mse=sum(se)/L
  time=(proc.time()-ptm)[1]
  #summarize
  #ape=sum(pe)/L
  if(model=="RA"||model=="RB")
  {
    pro=0
  }
  else 
  {
    pro=0.1
  }
  s=OutlierSummary(w,pro)
  #return
  list(pe=pe,time=time,w=w,beta=b,index1=index1,index2=index2,iter=iter,mse=mse,s=s)
}