simulate=function(L,n,beta=NULL,model=c("A","B","C","D"),p=NULL,method="PWLQ",matlab=NULL,seed=2014,useDataFile=FALSE,
                  standardize=FALSE,penalty1="1-w0",updateInitial=FALSE,criterion="BIC",intercept=FALSE,initial="plain",
                  range="cross")
{
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
  b=array(0,dim=c(L,p))
  w=array(0,dim=c(L,n))
  iter=rep(0,L)
  error=NULL
  sumofy=NULL
  #out=GenerateDataByModel(n=n,beta=beta,model=model)
  #evaluate(matlab,"rng(2015);")
  for(i in 1:L)
  {
    if(useDataFile)
    {
      out=data[[i]]
    }
    else
    {
      out=GenerateDataByModel(n=n,beta=beta,model=model)
    }
    
    #evaluate(matlab,"[X y]=GenerateData()")
    #xm=getVariable(matlab,"X")
    #ym=getVariable(matlab,"y")
    #out$x=xm$X
    #out$y=ym$y
    if(method=="LAD")
    {
      init=InitParam(out$x,out$y,method="LAD")
      b[i,]=init$beta
    }
    else if(method=="ROSS")
    {
      setVariable(matlab, X=out$x)
      setVariable(matlab, y=out$y)
      evaluate(matlab,"[betaRoss]=RossSimulate(X,y)")
      betaRoss=getVariable(matlab, "betaRoss")
      b[i,]=as.vector(betaRoss$betaRoss)
    }
    else if(method=="ADL")
    {
      init=InitParam(out$x,out$y,method="LAD")
      beta0=ifelse(init$beta==0,0.01,init$beta)
      w0=rep(0.99,n)
      res=rcdreg(out$x,out$y,penalty="ADL",nlambda1=2,nlambda2=100,beta0=beta0,w0=w0,delta=0.000001,maxIter=1000)
      index=BIC(res$wloss[1,],apply(matrix(res$beta[1,,],100,p)!=0+0,1,sum),n,p)
      
      b[i,]=res$beta[1,index,]
      #w[i,]=res$w
      #iter[i]=res$iter
    }
    else if(method=="LTS")
    {
      require(robustHD)
      res=sparseLTS(out$x,out$y,intercept=FALSE)
      b[i,]=res$coefficients
    }
    else if(method=="MMNNG")
    {
      #fbeta=paste("data\\",model,n,"X",p,"_beta.rda",sep="")
      #load(fbeta)
      #b=betaHat
      #break
      
      source("mmnngreg.R")
      res=mmnngreg(out$x,out$y)
      b[i,]=res$betac[-1]
    }
    else if(method=="MMNNG_DATA")
    {
      source("mmnngreg.R")
      res=mmnngreg(out$x,out$y)
      f=paste("data\\",model,n,"X",p,".rda",sep="")
      fbeta=paste("data\\",model,n,"X",p,"_beta.rda",sep="")
      lf=try(load(f))
      if(class(lf)=="try-error") # first time
      {
        data=list(out)
        betaHat=matrix(nrow=1,ncol=p,res$betac[-1])
      }
      else
      {
        if(length(data)==L)
          stop("data collect finish!")
        data=c(data,list(out))
        load(fbeta)
        betaHat=rbind(betaHat,res$betac[-1])
      }
      save(data,file=f)
      save(betaHat,file=fbeta)
    }
    else if(method=="COMPARE")
    {
      res1=mmnngreg(out$x,out$y)
      if(initial=="LTS")
      {
        require(robustHD)
        init=sparseLTS(out$x,out$y,intercept=intercept)
        beta0=SetBeta0(init$coefficients)
       
        w0=ifelse(init$wt==1,0.99,0.01)
      }
      else if(initial=="plain")
      {
        beta0=rep(1,p)
        w0=rep(0.99,n)
      }
      else
      {
        beta0=rep(1,p)
        w0=rep(0.99,n)
      }
      
      res=srcdreg(out$x,out$y,penalty1=penalty1,nlambda1=50,nlambda2=100,beta0=beta0,w0=w0,delta=0.000001,maxIter=1000,
                  intercept=intercept,standardize=standardize,updateInitial=updateInitial,criterion=criterion)
      return (list(beta=res$beta,w=res$w,beta1=res1$betac[-1]))
    }
    else
    {
      if(initial=="LTS")
      {
        require(robustHD)
        init=sparseLTS(out$x,out$y,intercept=intercept)
        beta0=SetBeta0(init$coefficients)
        w0=UpdateWeight(out$x,out$y,beta0)
        w0=ifelse(as.vector(w0)==1,0.99,w0)
      }
      else if(initial=="plain")
      {
        beta0=rep(1,p)
        w0=rep(0.99,n)
      }
      else
      {
        beta0=rep(1,p)
        w0=rep(0.99,n)
      }
      if(range=="cross")
      {
        res=srcdreg(out$x,out$y,penalty1=penalty1,nlambda1=50,nlambda2=100,beta0=beta0,w0=w0,delta=0.000001,maxIter=1000,
                    intercept=intercept,standardize=standardize,updateInitial=updateInitial,criterion=criterion)
      }
      else if(range=="all")
      {
        res=rcdreg(out$x,out$y,penalty1=penalty1,nlambda1=50,nlambda2=100,beta0=beta0,w0=w0,delta=0.000001,maxIter=1000,
                    intercept=intercept,standardize=standardize,updateInitial=updateInitial,criterion=criterion)
      }
                 
      
      #test rcdreg
      #res2=rcdreg(out$x,out$y,penalty="ADL",nlambda1=50,nlambda2=100,beta0=beta0,w0=w0,delta=0.000001,maxIter=1000)
     #final=BICPWLQ(res2$wloss,res2$beta,w=res2$w,res2$lambda1,res2$lambda2,n)
      b[i,]=res$beta
      w[i,]=res$w
      iter[i]=res$iter
      #b[i,]=final$beta
    }
  }
  
  #SE
  se=apply((t(b)-beta)^2,2,sum)
  
  #vs
  vs=rep(0,3)
  pNum=sum(beta!=0)
  bb=(abs(b)>=5e-5)
  #bb=(b>=5e-5)
  n1=apply(bb,1,sum)
  bb=matrix(bb[,1:pNum],nrow=L,ncol=pNum)
  n2=apply(bb,1,sum)
  vs[1]=sum((n1[n2==pNum]==pNum)+0)/L
  vs[2]=sum((n1[n2==pNum]>=pNum)+0)/L
  vs[3]=sum(n1)/L
  #return
  list(beta=b,se=se,vs=vs,iter=iter,w=w,sum=sumofy,error=error)
}