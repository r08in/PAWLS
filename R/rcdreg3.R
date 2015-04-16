#without search the whole grid of lambda
RCDReg3=function(x, y,penalty1=c("log","1-w0"),penalty2=c("MCP", "SCAD", "ADL"),lambda1,lambda2,beta0,w0,delta, 
                 maxIter,intercept=TRUE,updateInitial=TRUE)
{  
  #lambda2=lambda2[round(length(lambda2)/2):length(lambda2)]
  #lambda1=lambda1[1:length(lambda1)-1]
  
  L1=length(lambda1)
  L2=length(lambda2)
  m=dim(x)[2]
  n=length(y)
  pre1=pre2=0
  lmaxIter=10
  
  #initial lambda
  index2=L2/2
  res=InnerReg(x,y,penalty1=penalty1,penalty2=penalty2,lambda1,lambda2[index2],beta0,w0,delta,maxIter,intercept=intercept) #fix lambda2
  index1=BIC(as.vector(res$wloss),apply(matrix(res$w,L1,n)!=1+0,1,sum)+apply(matrix(res$beta,L1,m)!=0+0,1,sum),n,m,type="w") #find best lambda1
  #index1=BIC(as.vector(res$wloss),dfs(x,matrix(res$beta,L1,m),matrix(res$w,L1,n)),n,m,type="w") #find best lambda1
  res=InnerReg(x,y,penalty1=penalty1,penalty2=penalty2,lambda1[index1],lambda2,beta0,w0,delta,maxIter,intercept=intercept) #fix lambda1
  #index2=BIC(as.vector(res$wloss),dfs(x,matrix(res$beta,L2,m),matrix(res$w,L2,n)),n,m) #find best lambda2
  index2=BIC(as.vector(res$wloss),apply(matrix(res$w,L2,n)!=1+0,1,sum)+apply(matrix(res$beta,L2,m)!=0+0,1,sum),n,m) #find best lambda2
  
  ##loop to estimate and find the best
  iter=0
  
  while((pre1!=index1||pre2!=index2)&&(iter<lmaxIter))
  {
    iter=iter+1
    pre1=index1
    pre2=index2
    
    if(updateInitial)
    {
      beta0=matrix(res$beta,L2,m)[index2,]    
      w0=matrix(res$w,L2,n)[index2,]
      beta0=SetBeta0(beta0)
      w0=ifelse(w0==1,0.99,w0)
    }    
    res=InnerReg(x,y,penalty1=penalty1,penalty2=penalty2,lambda1,lambda2[index2],beta0,w0,delta,maxIter,intercept=intercept) #fix lambda2
    #index1=BIC(as.vector(res$wloss),dfs(x,matrix(res$beta,L1,m),matrix(res$w,L1,n)),n,m,type="w") #find best lambda1  
    index1=BIC(as.vector(res$wloss),apply(matrix(res$w,L1,n)!=1+0,1,sum)+apply(matrix(res$beta,L1,m)!=0+0,1,sum),n,m,type="w") #find best lambda1  
    
    if(updateInitial)
    {
      beta0=matrix(res$beta,L1,m)[index1,]
      w0=matrix(res$w,L1,n)[index1,]
      beta0=SetBeta0(beta0)
      w0=ifelse(w0==1,0.99,w0)
    }    
    res=InnerReg(x,y,penalty1=penalty1,penalty2=penalty2,lambda1[index1],lambda2,beta0,w0,delta,maxIter,intercept=intercept) #fix lambda1
    #index2=BIC(as.vector(res$wloss),dfs(x,matrix(res$beta,L2,m),matrix(res$w,L2,n)),n,m) #find best lambda2
    index2=BIC(as.vector(res$wloss),apply(matrix(res$w,L2,n)!=1+0,1,sum)+apply(matrix(res$beta,L2,m)!=0+0,1,sum),n,m) #find best lambda2
    if(pre2==index2&&pre1==index1) break;
  }
  
  if(updateInitial)
  {
    beta0=matrix(res$beta,L2,m)[index2,]    
    w0=matrix(res$w,L2,n)[index2,]
    beta0=SetBeta0(beta0)
    w0=ifelse(w0==1,0.99,w0)
  }  
  #test##
  res=InnerReg(x,y,penalty1=penalty1,penalty2=penalty2,lambda1[index1],lambda2[index2],beta0,w0,delta,maxIter,intercept=intercept) #fix lambda2 and lambda2
  #return 
  i=index1
  j=index2
  list(lambda1=lambda1[i],lambda2=lambda2[j],
       beta=as.vector(res$beta),w=as.vector(res$w),
       wloss=res$wloss,bdf=sum(res$beta!=0+0),wdf=sum(res$w!=1+0),
       index1=i,index2=j,iter=iter)
  
}

InnerReg=function(x, y,penalty1="1-w0",penalty2="ADL",lambda1,lambda2,beta0,w0,delta, maxIter,intercept=TRUE)
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
    loss[l1,1]=t(y)%*%y ##initial loss[l1,1]
    if(penalty1=="log")
    {
      lam1=sqrt((lambda1[l1]/abs(log(w0)))*n) #init sqrt(lambda1/abs(log(w0))n)
    }
    else #1-wo
    {
      lam1=(lambda1[l1]/abs(1-w0))*n
    }
    
    ##iteration for each lamda2
    for(l2 in lstart2:L2)
    {  
      lam2=lambda2[l2]/abs(beta0)
      if(intercept)
      {
        lam2[1]=0
      }
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
          ##(2)update betaj
          if (penalty2=="ADL")
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
        if(penalty1=="log")
        {
          absr=abs(r)
          w[l1,l2,]=ifelse(absr>lam1,lam1/absr,1)
        }
        else #1-w0
        {
          sqr=r^2
          w[l1,l2,]=ifelse(sqr>lam1,lam1/sqr,1)
        }
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

SetBeta0=function(beta0)
{
  if(sum(beta0)==0)
  {
    b=0.01
  }
  else
  {
    b=min(abs(beta0[beta0!=0]))*0.01
  }
  ifelse(beta0==0,b,beta0)
}
