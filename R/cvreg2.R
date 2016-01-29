cvreg=function(x, y,lambda1,lambda2,beta0,w0,K=5)
{
  
  L=2
  i=1
  L1=length(lambda1)
  L2=length(lambda2)
  n=length(y)
  p=dim(x)[2]
  ws=matrix(1,nrow=L1,ncol=n)
  bs=matrix(0,nrow=L2,ncol=p)
  size=floor(n/K)
  rseq=sample(1:n,n)
  
  while (i<=L) 
  {
    ##find lambda1 and update w by BIC
    #compute residule
    r=y-x%*%beta0
    sqr=r^2
    
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
      
      #compute pe
      for(l2 in 1:L2)
      {
        H=t(trainx)%*%trainx+n*lambda2[l2]*diag(1,p)
        bs[l2,]=solve(H)%*%t(trainx)%*%trainy
        r=testy-testx%*%bs[l2,]
        pe[l2]=pe[l2]+sum(r^2)/size
      }
      
    }
    pe=pe/K
    #find best lambda2 that have testPE smallset
    index2=which.min(pe)
    beta=bs[index2,]
    ##repeat 
    beta0=beta
    w0=w
    i=i+1
  }
  
  list(beta=beta,w=w)
  
}