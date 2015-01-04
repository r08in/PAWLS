
n=100
p=200
pNum=5
out=GenerateData (n,p,pNum,errorSigma=1,outlier.op="MEANSHIFT",outlier.pro=0.1,outlier.r=10)
final=BetaInit(out$x,out$y,method="LASSO")
