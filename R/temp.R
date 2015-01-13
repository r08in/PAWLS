GCDReg2=function(x, y, groupInfo, penalty,gamma, lamda, delta, maxIter)
{  
  ##declaration
  n=length(y)
  m=dim(x)[2]
  p=length(groupInfo)
  L=length(lamda)
  begin=rep(0,p)
  end=rep(0,p)
  lstart=1 ##since lamda[1] will give all beta 0
  
  ##reslut to be returned
  beta=matrix(ncol=m,nrow=L,0)
  loss=rep(0,L)
  iter=rep(0,L)
  
  ##temp
  r=rep(0,n)
  betaPre=rep(0,m)
  betaShift=rep(0,m)
  z=rep(0,m)
  tempb=rep(0,m)
  
  ##initial
  r=y    
  for(j in 1:p) ##initial z
  {
    begin=end+1
    end=end+groupInfo[j]
    tempz=CrossProduct2(x,y,begin,end)  ## division by n included
    for(i in begin:end)
    {
      z[i]=tempz[i-begin+1]
    }
  }
  
  loss[1]=t(y)%*%y ##initial loss[1]
  
  ##iteration for each lamda
  for(l in lstart:L)
  {  
    if(l>=2)
    {
      betaPre=beta[l-1,] ##assign previous beta to betaPre
    }
    
    ##iteration for all covariates
    while(iter[l]<maxIter)
    {
      iter[l]=iter[l]+1
      begin=end=0
      
      ##iteration for each covariate group
      for(j in 1:p)
      {
        begin=end+1
        end=end+groupInfo[j]
        
        ##(1)calculate z j group
        tempz=CrossProduct2(x,r,begin,end)
        for(i in begin:end)
        {
          z[i]=tempz[i-begin+1]+betaPre[i]
        }
        
        ##(2)update beta j group
        if (penalty=="MCP")
        {
          tempb=McPGroup2(z,begin,end,lamda[l],gamma)
          beta[l,begin:end]=tempb
        }
        
        ##(3)update r
        betaShift[begin:end]=beta[l,begin:end]-betaPre[begin:end]
        r=r-x[,begin:end]%*%betaShift[begin:end]  
        
      }
      
      ##update betaPre for next iteration
      betaPre=beta[l,]
      
      ## Check for convergence
      if(t(betaShift)%*%betaShift<delta)
      {
        break;
      }
      
    }
    ##compute square of loss
    loss[l]=t(r)%*%r
  }
  
  list(beta,loss,iter)
}