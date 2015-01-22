
n=100
p=200
pNum=8
oPro=0.1
beta=c(3,1.5,0,0,2,0,0,0)
#-----------------model1-----------------------------------####
#out=GenerateData(n=n,dataSetNum=1,beta=beta,errorSigma=2,r=0.5)
#-----------------model2-----------------------------------####
out=GenerateData(n=n,dataSetNum=1,beta=beta,errorSigma=2,r=0.95)


#out=GenerateData (n,p,pNum,errorSigma=1,outlier.op="MEANSHIFT",outlier.pro=oPro,outlier.r=10)
init=InitParam(out$x,out$y,method="LAD")
beta0=ifelse(init$beta==0,0.01,init$beta)
w0=ifelse(init$weight==1,0.99,init$weight)
res=rcdreg(out$x,out$y,penalty="ADL",nlambda1=50,nlambda2=100,
      beta0=beta0,w0=w0,delta=0.000001,maxIter=1000)


final=BICPWLQ2(res$wloss,res$beta,w=res$w,res$lambda1,res$lambda2,n)
s=MySummary(res$beta,res$w,pNum,round(n*oPro))
b=s$b
ww=s$w

bic=final$res$bic
save(res,file = "res.rda")

###########
outcome=NULL
for(i in 1:10)
{
  n=100
  p=200
  pNum=10
  out=GenerateData (n,p,pNum,errorSigma=1,outlier.op="MEANSHIFT",outlier.pro=0.1,outlier.r=10)
  init=InitParam(out$x,out$y,method="LAD")
  beta0=ifelse(init$beta==0,sum(init$beta)/p,init$beta)
  w0=ifelse(init$weight==0,sum(init$weight)/n,init$weight)
  res=rcdreg(out$x,out$y,penalty="ADL",nlambda1=50,nlambda2=200,
             beta0=beta0,w0=w0,delta=0.000001,maxIter=1000)
  final=BICPWLQ(res$wloss,res$beta,res$w,res$lambda1,res$lambda2,n)
  outcome=c(outcome,list(final))
}

#####integrative study#####
n=100
p=200
pNum=5
dataSetNum=4
nlambda=100
groupInfo=rep(dataSetNum,p)
data=GenerateData(n=n,p=p,pNum=pNum,dataSetNum=dataSetNum,errorSigma=1,
                  offSet=0,outlier.op="MEANSHIFT",outlier.pro=0.1,outlier.r=10)#outlier.op="MEANSHIFT"
out=CombineMultiLm(data$x,data$y)
init=InitParam(out$x,out$y,method="LAD")
w0=ifelse(init$weight==1,0.99,init$weight)
#out=data
res=rgcdreg(out$x,out$y,groupInfo,penalty="MCP",gamma=3,nlambda1=50,nlambda2=100,
       w0=w0,delta=0.000001,maxIter=1000)
final=BICPWLQ(res$wloss,res$beta,res$w,res$lambda1,res$lambda2,n*dataSetNum)
s=MySummary(res$beta,res$w,pNum*dataSetNum,round(n*oPro))
b=s$b
w=s$w

#####model 1#####
beta=c(3,1.5,0,0,2,0,0,0)
GenerateData(n=n,p=p,pNum=pNum,dataSetNum=1,beta=beta,errorSigma=2,r=0.5)