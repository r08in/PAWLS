

MySummary=function(beta,weight,pnum,onum)
{
  l1=dim(beta)[1]
  l2=dim(beta)[2]
  b=matrix(0,nrow=l1,ncol=l2)
  w=matrix(0,nrow=l1,ncol=l2)
  for(i in 1:l1)
    for(j in 1:l2)
    {
      b[i,j]=sum(beta[i,j,]!=0+0)+sum(beta[i,j,(1:pnum)]!=0+0)*0.01
      w[i,j]=sum(weight[i,j,]!=1+0)+sum(weight[i,j,(1:onum)]!=1+0)*0.01
    }
  
  list(b=b,w=w)
}

MySummary2=function(beta,weight,pnum,onum)
{
  l1=length(beta)
  l2=dim(beta[[1]])[1]
  b=matrix(0,nrow=l1,ncol=l2)
  w=matrix(0,nrow=l1,ncol=l2)
  for(i in 1:l1)
    for(j in 1:l2)
    {
      b[i,j]=sum(beta[[i]][j,]!=0+0)+sum(beta[[i]][j,(1:pnum)]!=0+0)*0.01
      w[i,j]=sum(weight[[i]][j,]!=1+0)+sum(weight[[i]][j,1:onum]!=1+0)*0.01
    }
  
  list(b=b,w=w)
}