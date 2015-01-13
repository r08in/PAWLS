## This function is to combine linear model from each dataset
## To form a group linear model

CombineMultiLm=function(x,y)
{
  #error checking
  if(length(dim(x))!=3)
    stop("The dimension of multi-design matrix x should be 3.")
  if(length(dim(y))!=2)
    stop("The dimension of multi-observations y should be 2.")
  if(dim(x)[1]<=1)
    stop("The number of dataset should greater than 1.")
  if(dim(x)[1]!=dim(y)[1]||dim(x)[2]!=dim(y)[2])
    stop("The dimensions of x and y are not matched. ")
  
  m=dim(x)[1]
  n=dim(x)[2]
  p=dim(x)[3]
  
  #combine Y
  tempy=NULL
  for(i in 1:m)
  {
    tempy=c(tempy,y[i,])
  }  
  
  #combine x
  tempx=NULL  
  for(j in 1:p)
  {
    xx=NULL
    for(i in 1:m)
    {
      xx=cbind(xx,c(rep(0,(i-1)*n),x[i,,j],rep(0,(m-i)*n)))
    }
    tempx=cbind(tempx,xx)
  }
  
  list(x=tempx,y=tempy)
}


##This function is to combine multiple dataset
##each dataset is already in a group linear model

CombineDataset=function(XX,xinfo,pinfo,YY)
{
  ##checking data
  
  ##combine x
  m=length(xinfo)
  p=length(pinfo)
  colStart=1
  X=NULL
  for(j in 1:p) #for each Xj
  {    
    colSum=pinfo[j]*m
    Xj=matrix(0,0,colSum)
    leftColNum=0
    rowStart=1
    
    for(i in 1:m) #for combine each xj in dataset i
    {
      tempx=cbind(matrix(0,xinfo[i],leftColNum),
                  XX[rowStart:(rowStart+xinfo[i]-1),colStart:(colStart+pinfo[j]-1)],
                  matrix(0,xinfo[i],colSum-leftColNum-pinfo[j]))
      Xj=rbind(Xj,tempx)
      leftColNum=leftColNum+pinfo[j]
      rowStart=rowStart+xinfo[i]
    }
    colStart=colStart+pinfo[j]
    X=cbind(X,Xj)
  }
  
  ##combine y(assume YY has been combine)
  
  ##return list
  list(x=X,y=YY)
}


