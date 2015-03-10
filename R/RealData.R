#This file for real data analysis

PrepareData=function(file="fltr.csv",yname="1389163_at",screenNum=500)
{
  data=read.table(file,header=TRUE)
  data=t(data)
  y=data[-1,(data[1,]==yname)]
  x=data[-1,(data[1,]!=yname)]
  #colnames(x)=data[1,(data[1,]!=yname)]
  y=as.numeric(y)
  x=apply(x,2,as.numeric)
  x=ScreenData(y,x,screenNum)$x
  list(y=y,x=x,data=data)
  
}

ScreenData=function(y,x,screenNum)
{
  m=dim(x)[2]
  snum1=ifelse(m>4*screenNum,4*screenNum,m)
  index1=order(apply(x,2,var),decreasing=TRUE)[1:snum1]
  xx=x[,index1]
  index=order(abs(cor(y,xx)),decreasing=TRUE)[1:screenNum]
  list(y=y,x=xx[,index])
}

AddIntercept=function(x)
{
  n=dim(x)[1]
  x1=rep(1,n)
  xx=cbind(x1,x)
  xx
}