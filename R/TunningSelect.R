## This file is for tunning parameter selection

##BIC

BICSelect2=function(rss,n,betas,lambdas,inv=1)
{
  #declare and initial
  m=length(lambdas)
  indexSelected=1
  df=sum(betas[1,]!=0)
  BicPre=log(rss[1]/n)+log(n)*df/n
  res=matrix(0,m,3)
  res[1,1]=lambdas[1]
  res[1,2]=BicPre
  res[1,3]=df
  BicPre=1000000 # not use the  first lambda(all zero)
  start=round(m*((1-inv)/2))
  end=round(m*(0.5+inv/2))
  for(i in 2:m)
  {
    df=sum(betas[i,]!=0)
    BicTemp=log(rss[i]/n)+log(n)*df/n
    res[i,1]=lambdas[i]
    res[i,2]=BicTemp
    res[i,3]=df
    if(i>=start&&i<=end&&BicTemp<=BicPre)
    {
      indexSelected=i
      BicPre=BicTemp
    }
  }
  res=data.frame(lambda=res[,1],BIC=res[,2],df=res[,3])
  list(lambda=lambdas[indexSelected],beta=betas[indexSelected,],rss=rss[indexSelected],
       df=sum(betas[indexSelected,]!=0),index=indexSelected,res=res,betas=betas)
}

GreedySelect=function(betas,lambdas,inv=0.9)
{
  m=length(lambdas)
  index=round(m*(0.5+inv/2))
  list(lambda=lambdas[index],df=sum(betas[index,]!=0),index=index,beta=betas[index,])
}

BICPWLQ=function(wloss,beta,w,lambda1,lambda2,n,inv=1)
{
  #declare and initial
  l1=length(lambda1)
  l2=length(lambda2)
  index1=index2=1
  bicPre=100000000
  bicTemp=matrix(0,l1,l2)
  wdf=matrix(0,l1,l2)
  bdf=matrix(0,l1,l2)
  start1=round(l1*((1-inv)/2))+1
  end1=round(l1*(0.5+inv/2))
  start2=round(l2*((1-inv)/2))+1
  end2=round(l2*(0.5+inv/2))
  for(i in start1:end1)
  {
    for(j in start2:end2)
  {
    wdf[i,j]=sum(w[i,j,]!=1+0)
    bdf[i,j]=sum(beta[i,j,]!=0+0)
    bicTemp[i,j]=log(wloss[i,j]/(n)+1)+(bdf[i,j]+wdf[i,j])*log(n)/(n)
    if(bicTemp[i,j]<=bicPre)
    {
      index1=i
      index2=j
      bicPre=bicTemp[i,j]
    }
  }
  }
  res=list(lambda1=lambda1,lambda2=lambda2,bdf=bdf,wdf=wdf,bic=bicTemp)
  i=index1
  j=index2
  list(lambda1=lambda1[i],lambda2=lambda2[j],
       beta=beta[i,j,],w=w[i,j,],
       wloss=wloss[i,j],bdf=bdf[i,j],wdf=wdf[i,j],
       index1=i,index2=j,res=res)
}

BICPWLQ2=function(wloss,beta,w,lambda1,lambda2,n,inv=1)
{
  #declare and initial
  l1=length(lambda1)
  l2=length(lambda2)
  p=dim(beta)[3]
  index1=index2=1
  bicPre=100000000
  bicTemp=matrix(0,l1,l2)
  wdf=matrix(0,l1,l2)
  bdf=matrix(0,l1,l2)
  start1=round(l1*((1-inv)/2))+1
  end1=round(l1*(0.5+inv/2))
  start2=round(l2*((1-inv)/2))+1
  end2=round(l2*(0.5+inv/2))
  index1=l1-1
  index2=BIC(wloss[index1,],apply(beta[index1,,]!=0+0,1,sum),n,p,type="beta")
  index1=BIC(wloss[,index2],apply(w[,index2,]!=1+0,1,sum),n,p,type="w")
  index2=BIC(wloss[index1,],apply(beta[index1,,]!=0+0,1,sum),n,p,type="beta")
  res=list(lambda1=lambda1,lambda2=lambda2,bdf=bdf,wdf=wdf,bic=bicTemp)
  i=index1
  j=index2
  list(lambda1=lambda1[i],lambda2=lambda2[j],
       beta=beta[i,j,],w=w[i,j,],
       wloss=wloss[i,j],bdf=sum(beta[index1,index2,]!=0+0),wdf=sum(w[index1,index2,]!=1+0),
       index1=i,index2=j,res=res)
}

BIC=function(loss,df,n,p,type="beta")
{
  if(type=="beta")
  {
    vl=(log(loss/n+1)+log(n)*df/(n))
    #vl=loss/n+log(n)*df/(n)
    #vl=loss/n+log(n)*df/(n)+log(choose(p,df))/n
  }
  else 
  {
    vl=(log(loss/n+1)+log(n)*df/(n))
    #vl=loss/n+log(n)*df/(n)
    #vl=loss/n+log(n)*df/(n)+log(choose(p,df))/n
  }

 
  l=length(df)
  index=(1:l)[vl==min(vl)]
  if(length(index)!=1)
  {
    index[1]
  }
  else
  {
    index
  }
}

dfs=function(x,beta,w)
{
  if(dim(beta)[1]!=dim(w)[1])
    stop("in df: beta and w should have the same row!")
  require('Matrix')
  m=dim(beta)[1]
  df=rep(0,m)
  b=abs(beta)>5e-5
  for(i in 1:m)
  {
    if(sum(b[i,])==0) next
    df[i]=rankMatrix(as.matrix(w[i,]*x[,b[i,]]))
  }
  df
}