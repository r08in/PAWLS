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