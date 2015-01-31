
n=50
p=8
pNum=8
beta1=c(3,2,1.5,1,1,1,1,1)
beta2=c(3,2,1.5,0,0,0,0,0)
#-----------------model1(normal)-----------------------------------####
#out=GenerateData(n=n,dataSetNum=1,beta=beta,errorSigma=2,r=0.5)
#-----------------model2(heavy tail)-----------------------------------####
 out=GenerateData(n=n,dataSetNum=1,beta=beta,errorSigma=2,errorType="t",r=0.5)


out=GenerateData (n,p,pNum=1,errorSigma=1,outlier.op="MEANSHIFT",outlier.pro=oPro,outlier.r=10)
init=InitParam(out$x,out$y,method="LAD")
beta0=ifelse(init$beta==0,0.01,init$beta)
w0=ifelse(init$weight==1,0.99,init$weight)
res=rcdreg(out$x,out$y,penalty="ADL",nlambda1=50,nlambda2=100,
      beta0=beta0,w0=w0,delta=0.000001,maxIter=1000)


final=BICPWLQ2(res$wloss,res$beta,w=res$w,res$lambda1,res$lambda2,n)
s=MySummary(res2$beta,res2$w,3,round(n*0.1))
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

#simulation for 4 different linear model

L=100
n=50
p=8
beta1=c(3,2,1.5,1,1,1,1,1)
beta2=c(3,2,1.5,0,0,0,0,0)
beta=beta2
#out=GenerateDataByModel(n=n,beta=beta,model="A")
#out=GenerateDataByModel(n=n,beta=beta,model="B")
#out=GenerateDataByModel(n=n,beta=beta,model="C")
out=GenerateDataByModel(n=n,beta=beta,model="D")
init=InitParam(out$x,out$y,method="LAD")
beta0=ifelse(init$beta==0,0.01,init$beta)
w0=ifelse(init$weight==1,0.99,init$weight)
res=rcdreg(out$x,out$y,penalty="ADL",nlambda1=50,nlambda2=100,
           beta0=beta0,w0=w0,delta=0.000001,maxIter=1000)

#---------------------simulation for 4 moedl---------------------------#
#--------------------compute se,vs(cs,cr,an)---------------------------#
#matLabDir="C:\\Users\\Administrator\\Desktop\\Sim"
matLabDir="D:\\matlab\\Sim"
matlab=PrepareMatlab(matLabDir)
#setting
L=50
n=50
p=8
beta1=c(3,2,1.5,1,1,1,1,1)
beta2=c(3,2,1.5,0,0,0,0,0)
beta=beta2
#simulation
out1=simulate(L,n,beta,"A")
out2=simulate(L,n,beta,"B")
out3=simulate(L,n,beta,"C")
out4=simulate(L,n,beta,"D")
outRoss4=simulate(L,n,beta,model="D",method="ROSS",matlab=matlab)
outRoss3=simulate(L,n,beta,model="C",method="ROSS",matlab=matlab)
outRoss2=simulate(L,n,beta,model="B",method="ROSS",matlab=matlab)
outRoss1=simulate(L,n,beta,model="A",method="ROSS",matlab=matlab)

#------------compare square error by boxplot----------------------#
se1=c(out1$se,outRoss1$se)
se2=c(out2$se,outRoss2$se)
se3=c(out3$se,outRoss3$se)
se4=c(out4$se,outRoss4$se)
group=c(rep("PWLQ",L),rep("ROSS",L))
m1=data.frame(se=se1,group=group)
m2=data.frame(se=se2,group=group)
m3=data.frame(se=se3,group=group)
m4=data.frame(se=se4,group=group)
boxplot(se~group,data=m1,main="model A")
x11();
boxplot(se~group,data=m2,main="model B")
x11();
boxplot(se~group,data=m3,main="model C")
x11();
boxplot(se~group,data=m4,main="model D")



#----------------another setting---------------------#

#setting
L=100
n=100
p=500
beta=c(rep(1,10),rep(0,p-10))

#simulation
out1=simulate(L,n,beta,"A")
out2=simulate(L,n,beta,"B")
out3=simulate(L,n,beta,"C")
out4=simulate(L,n,beta,"D")

out1_500=out1
out2_500=out2
out3_500=out3
out4_500=out4

