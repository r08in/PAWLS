
n=100
p=200
pNum=5
out=GenerateData (n,p,pNum,errorSigma=1,outlier.op="MEANSHIFT",outlier.pro=0.1,outlier.r=10)
init=InitParam(out$x,out$y)
beta0=ifelse(init$beta==0,0.001,init$beta)
w0=ifelse(init$weight==0,0.001,init$weight)
rcdreg(out$x,out$y,penalty="ADL",nlambda1=100,nlambda2=100,
      beta0=beta0,w0=w0,delta=0.000001,maxIter=1000)
init$beta
