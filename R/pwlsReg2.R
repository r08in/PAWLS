pwlsReg2=function(x, y,lambda,w0,delta, maxIter,intercept=TRUE)
{
  
  ##declaration
  n=length(y)
  m=dim(x)[2]
  p=m
  L=length(lambda)
  lstart=1 ##should be 1
  
  ##reslut to be returned
  beta=matrix(nrow=L,ncol=m,0)
  w=matrix(nrow=L,ncol=n,1)
  loss=rep(0,L)
  wloss=rep(0,L)
  iter=rep(0,L)
  
  betaPre=rep(0,m)
  wPre=rep(1,n)
  r=y
  ##iteration for each lamda1
  for(l1 in lstart:L)
  {
    ##initial
    
    shift=rep(0,n)
    lam1=(lambda[l1]/abs(1-w0))*n
    ##iteration for each lamda2
      ##iteration for all covariates
      while(iter[l1]<maxIter)
      {
        iter[l1]=iter[l1]+1
        
        #calculate coefficient c
        c=apply((x*wPre)^2,2,sum)/n
        
        ##iteration for each beta
        for(j in 1:m)
        {        
          ##(1)calculate zj 
          zj=t(x[,j]*wPre^2)%*%r/n+c[j]*betaPre[j]
          
          ##(2)update beta
          beta[l1,j]=UpdateBeta(zj,0,c[j])
          
          ##(3)update r
          r=r-x[,j]*(beta[l1,j]-betaPre[j])  
        }
        
        ##update w
        sqr=r^2
        w[l1,]=ifelse(sqr>lam1,lam1/sqr,1)
        shift=w[l1,]-wPre
        
        ##update betaPre and wPre for next iteration
        betaPre=beta[l1,]
        wPre=w[l1,]
        ## Check for convergence
        if(max(abs(shift))<delta)
        {
          break;
        }
        
      } #end for the inner loop
      
      ##compute square of loss
      loss[l1]=t(r)%*%r
      wloss[l1]=t(r*wPre)%*%(r*wPre)
    }#end iteration for  lambda
  
  list(beta=beta,loss=loss,iter=iter,w=w,wloss=wloss)
}
