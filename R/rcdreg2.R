RCDReg2=function(x, y,penalty,lambda1,lambda2,beta0,w0,delta, maxIter)
{  
  ##declaration
  n=length(y)
  m=dim(x)[2]
  p=m
  L1=length(lambda1)
  L2=length(lambda2)
  lstart1=1 ##should be 1
  lstart2=1
  
  ##reslut to be returned
  beta=array(0,dim=c(L1,L2,m))
  w=array(1,dim=c(L1,L2,n))
  loss=matrix(ncol=L2,nrow=L1,0)
  wloss=matrix(ncol=L2,nrow=L1,0)
  iter=matrix(ncol=L2,nrow=L1,0)
  
  betaPre=rep(0,m)
  wPre=rep(1,n)
  r=y
  ##iteration for each lamda1
  for(l1 in lstart1:L1)
  {
    ##initial
    
    shift=rep(0,m+n)
    c=rep(0,m)
    #z=t(x)%*%r/n#init z z=r%*%t(apply(x*wPre^2,2,sum)/n)
    loss[l1,1]=t(y)%*%y ##initial loss[l1,1]
    lam1=sqrt((lambda1[l1]/abs(log(w0)))*n) #init sqrt(lambda1/abs(log(w0))n)
    
    ##iteration for each lamda2
    for(l2 in lstart2:L2)
    {  
      lam2=lambda2[l2]/abs(beta0)
      
      ##iteration for all covariates
      while(iter[l1,l2]<maxIter)
      {
        iter[l1,l2]=iter[l1,l2]+1
        
        #calculate coefficient c
        c=apply((x*wPre)^2,2,sum)/n
        
        ##iteration for each beta
        for(j in 1:m)
        {        
          ##(1)calculate zj 
          zj=t(x[,j]*wPre^2)%*%r/n+c[j]*betaPre[j]
          if(is.na(zj))
          {
            cd=1
          }
          
          ##(2)update betaj
          if (penalty=="ADL")
          {
            beta[l1,l2,j]=UpdateBeta(zj,lam2[j],c[j])
            #beta[[l1]][l2,j]=UpdateBeta(zj,lam2[j]*sqrt(c[j]),c[j])
            #beta[[l1]][l2,j]=zj
          }
          
          ##(3)update r
          shift[j]=beta[l1,l2,j]-betaPre[j]
          r=r-x[,j]*shift[j]   
        }
        
        ##update w
        absr=abs(r)
        #test W=1
        # w[l1,l2,]=ifelse(absr>lam1,lam1/absr,1)
        w[l1,l2,]=rep(1,n)
        shift[(m+1):(m+n)]=w[l1,l2,]-wPre
        
        ##update betaPre and wPre for next iteration
        betaPre=beta[l1,l2,]
        wPre=w[l1,l2,]
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

UpdateBeta=function(zj,lambda,c=1)
{
  if(zj>lambda)
  {
    return((zj-lambda)/c)
  }
  else if(zj+lambda<0)
  {
    return((zj+lambda)/c)
  }
  else
  {
    return(0)
  }
}