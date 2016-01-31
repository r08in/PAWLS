cvreg2=function(x, y,lambda1=NULL,nlambda1=30,lambda2=NULL,nlambda2=30,beta0,w0,K=5,matlab=NULL)
{
  
  L=2
  i=1
  if(!is.null(lambda1))
  {
    nlambda1=length(lambda1)
  }
  L1=nlambda1
  n=length(y)
  p=dim(x)[2]
  ws=matrix(1,nrow=L1,ncol=n)
  size=floor(n/K)
  rseq=sample(1:n,n)
  
  while (i<=L) 
  {
    ##find lambda1 and update w by BIC
    #compute residule
    r=y-x%*%beta0[-1]-beta0[1]
    sqr=r^2
    
    #set up lambda1 for w
    temp=abs(1-w0)*sqr/n
    maxLambda1=max(temp)
    minLambda1=quantile(temp,probs=0.5)
    lambda1=seq(from=minLambda1,to=maxLambda1,length.out=L1)
    
    #update w for each lambda1
    for(l1 in 1:L1)
    {
      lam1=(lambda1[l1]/abs(1-w0))*n
      ws[l1,]=ifelse(sqr>lam1,lam1/sqr,1)
    }
    
    #find optimal lambda1 by BIC
    wloss=ws%*%sqr
    dfw=apply(ws!=1,1,sum)
    index1=BIC(loss=wloss,dfw=dfw,dfb=rep(0,L1),n=n,p=p,type="w")
    
    ##w
    w=ws[index1,]
    wy=w*y
    wx=diag(w)%*%x
    
    
    ##find lambda2 and update beta by CV
    if(is.null(lambda2))
    {
      lambda2=GetRidgeLambda(wx,wy,matlab=matlab)
    }
    L2=length(lambda2)
    bs=matrix(0,nrow=L2,ncol=p+1)
    #interation for each fold
    pe=rep(0,L2)
    for(k in 1:K)
    {
      #prepare data
      range=rseq[((k-1)*size+1):(k*size)]
      trainx=wx[-range,]
      trainy=wy[-range]
      testx=wx[range,]
      testy=wy[range]
      yy=c(trainy, rep(0,p))
      #compute pe
      for(l2 in 1:L2)
      {
        xx=cbind( c(rep(1,n-size),rep(0,p)), rbind(trainx,diag(sqrt(lambda2[l2]),p)))
        bs[l2,]=solve(t(xx)%*%xx)%*%t(xx)%*%yy
#         H=t(trainx)%*%trainx+n*lambda2[l2]*diag(1,p)
#         bs[l2,]=solve(H)%*%t(trainx)%*%trainy
        r=testy-testx%*%bs[l2,-1]-bs[l2,1]
        pe[l2]=pe[l2]+sum(r^2)/size
      }
      
    }
    pe=pe/K
    #find best lambda2 that have testPE smallset
    index2=which.min(pe)
    beta=bs[index2,]
    ##repeat 
    beta0=beta
    w0=ifelse(w==1,0.99,w)
    i=i+1
  }
  
  list(beta=beta,w=w,index1=index1,index2=index2,lambda1=lambda1,lambda2=lambda2)
  
}