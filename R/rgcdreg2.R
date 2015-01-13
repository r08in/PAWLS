RGCDReg2=function(x,y,groupInfo,penalty,gamma,lambda1,lambda2,w0,delta, maxIter)
{  
  ##declaration
  n=length(y)
  m=dim(x)[2]
  p=length(groupInfo)
  begin=rep(0,p)
  end=rep(0,p)
  L1=length(lambda1)
  L2=length(lambda2)
  lstart1=1 ##should be 1
  lstart2=1
  
  ##reslut to be returned
  beta=NULL
  w=NULL
  loss=matrix(ncol=L2,nrow=L1,0)
  wloss=matrix(ncol=L2,nrow=L1,0)
  iter=matrix(ncol=L2,nrow=L1,0)
  
  #initial
  betaPre=rep(0,m)
  wPre=rep(1,n)
  r=y
  for(i in 1:p) end[i]=sum(groupInfo[1:i])
  begin=c(1,end[1:p-1]+1)
  
  ##iteration for each lambda1
  for(l1 in lstart1:L1)
  {
    ##initial
    beta=c(beta,list(matrix(ncol=m,nrow=L2,0)))
    w=c(w,list(matrix(ncol=n,nrow=L2,1)))
    
    shift=rep(0,m+n)
    #z=t(x)%*%r/n#init z z=r%*%t(apply(x*wPre^2,2,sum)/n)
    #loss[l1,1]=t(y)%*%y ##initial loss[l1,1]
    lam1=sqrt((lambda1[l1]/abs(log(w0)))*n) #init sqrt(lambda1/abs(log(w0))n)
    
    ##iteration for each lambda2
    for(l2 in lstart2:L2)
    {  
      
      ##iteration for all covariates
      while(iter[l1,l2]<maxIter)
      {
        iter[l1,l2]=iter[l1,l2]+1
        xw=x*wPre
        ##iteration for each group beta        
        for(j in 1:p)
        {      
          jj=begin[j]:end[j]
          #calculate matrix coefficient cj
          cj=t(xw[,jj])%*%xw[,jj]/n #cj=xj'WWxj/n          
          
          ##(1)calculate zj 
          zj=t(x[,jj]*wPre^2)%*%r/n+cj%*%betaPre[jj]
          
          ##(2)update betaj
          if (penalty=="MCP")
          {
            beta[[l1]][l2,jj]=McPGroup(zj,groupInfo[j],lambda2[l2],gamma,cj)
          }
          
          ##(3)update r
          shift[jj]=beta[[l1]][l2,jj]-betaPre[jj]
          r=r-x[,jj]%*%shift[jj]   
        }
        
        ##update w
        absr=abs(r)
        w[[l1]][l2,]=ifelse(absr>lam1,lam1/absr,1)
        shift[m+1:m+n]=w[[l1]][l2,]-wPre
        
        ##update betaPre and wPre for next iteration
        betaPre=beta[[l1]][l2,]
        wPre=w[[l1]][l2,]
        
        ## Check for convergence
        if(t(shift)%*%shift<delta)
        {
          break;
        }
        
      } #end for the inner loop
      
      ##compute square of loss
      loss[l1,l2]=t(r)%*%r
      wloss[l1,l2]=t(r*wPre)%*%(r*wPre)
    }#end iteration for each lambda2 fixed lambda1
    
  }#end iteration for each lambda1
  
  list(beta=beta,loss=loss,iter=iter,w=w,wloss=wloss)
}


McPGroup=function(zj,groupSize,lambda,gamma,cj)
{
  ##error checking
  #if(end<begin||end<0||begin<0)
  invcj=ifelse(cj!=0,1/cj,0)
  sinvcj=sqrt(invcj)
  zzj=sinvcj%*%zj
  znorm=sqrt(t(zzj)%*%zzj)  
  lambda2=lambda*sqrt(groupSize)
  if(is.na(znorm))
  {
    a=1
  }
  if(znorm<=lambda2) ## set all to zero
  {
    tempb=rep(0,groupSize)
  }
  else if(lambda2<znorm && lambda2*gamma>=znorm)
  {
    tempb=(invcj%*%zj)%*%(gamma/(gamma-1)*(1-lambda2/znorm))
  }
  else
  {
    tempb=invcj%*%zj
  }
  tempb
}
