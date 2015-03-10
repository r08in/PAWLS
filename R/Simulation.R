simulate=function(L,n,beta,model=c("A","B","C","D"),method="PWLQ",matlab=NULL,seed=2014)
{
  if(!is.null(seed))
  {
    set.seed(seed)
  }
  p=length(beta)
  b=array(0,dim=c(L,p))
  w=array(0,dim=c(L,n))
  iter=rep(0,L)
  #out=GenerateDataByModel(n=n,beta=beta,model=model)
  #evaluate(matlab,"rng(2015);")
  for(i in 1:L)
  {
    out=GenerateDataByModel(n=n,beta=beta,model=model)
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
    else
    {
      #beta0=ifelse(init$beta==0,0.01,init$beta)
      #w0=ifelse(init$weight==1,0.99,init$weight)
      beta0=rep(1,p)
      w0=rep(0.99,n)
      res=srcdreg(out$x,out$y,penalty="ADL",nlambda1=50,nlambda2=100,beta0=beta0,w0=w0,delta=0.000001,maxIter=1000,intercept=FALSE)
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
  n2=apply(bb[,1:pNum],1,sum)
  vs[1]=sum((n1[n2==pNum]==pNum)+0)/L
  vs[2]=sum((n1[n2==pNum]>=pNum)+0)/L
  vs[3]=sum(n1)/L
  #return
  list(beta=b,se=se,vs=vs,iter=iter,w=w)
}