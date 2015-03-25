
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

L=2
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
L=100
n=50
p=8
beta1=c(3,2,1.5,1,1,1,1,1)
beta2=c(3,2,1.5,0,0,0,0,0)
beta=beta2
#simulation for pwls
out1=simulate(L,n,beta,"A",method="MMNNG")
SaveResult(out1$vs,"out1.txt")
out2=simulate(L,n,beta,"B")
SaveResult(out2$vs,"out2.txt")
out3=simulate(L,n,beta,"C",method="MMNNG")
SaveResult(out3$vs,"out3.txt")
out4=simulate(L,n,beta,"D",method="MMNNG",seed=20)
out4=simulate(L,n,beta,"D")
SaveResult(out4$vs,"out4.txt")

#simulation for SROS
matLabDir="D:\\matlab\\Sim"
matlab=PrepareMatlab(matLabDir)
outRoss4=simulate(L,n,beta,model="D",method="ROSS",matlab=matlab)
outRoss3=simulate(L,n,beta,model="C",method="ROSS",matlab=matlab)
outRoss2=simulate(L,n,beta,model="B",method="ROSS",matlab=matlab)
outRoss1=simulate(L,n,beta,model="A",method="ROSS",matlab=matlab)

#simulation for ADL

outADL1=simulate(L,n,beta,"A",method="ADL")
outADL2=simulate(L,n,beta,"B",method="ADL")
outADL3=simulate(L,n,beta,"C",method="ADL")
outADL4=simulate(L,n,beta,"D",method="ADL")

OulierSummary(out3$w)

out5=simulate(L,n,beta,"E")

outADL1=simulate(L,n,beta,"A",method="ADL")
outADL2=simulate(L,n,beta,"B",method="ADL")
outADL3=simulate(L,n,beta,"C",method="ADL")
outADL4=simulate(L,n,beta,"D",method="ADL")

#------------compare square error by boxplot----------------------#
se1=c(out1$se,outRoss1$se,outADL1$se)
se2=c(out2$se,outRoss2$se,outADL2$se)
se3=c(out3$se,outRoss3$se,outADL3$se)
se4=c(out4$se,outRoss4$se,outADL4$se)

se1=c(out1$se,outRoss1$se)
se2=c(out2$se,outRoss2$se)
se3=c(out3$se,outRoss3$se)
se4=c(out4$se,outRoss4$se)
group=c(rep("PWLQ",L),rep("SROS",L))
m1=data.frame(se=se1,group=group)
m2=data.frame(se=se2,group=group)
m3=data.frame(se=se3,group=group)
m4=data.frame(se=se4,group=group)
attach(mtcars)
par(mfrow=c(2,2))
boxplot(se~group,data=m1,main="model A")
boxplot(se~group,data=m2,main="model B")
boxplot(se~group,data=m3,main="model C")
boxplot(se~group,data=m4,main="model D",ylim=c(0,8))



#----------------another setting---------------------#

#setting
L=100
n=100
p=500
beta=c(rep(2,10),rep(0,p-10))
#simulation
out1=simulate(L,n,beta,"A")
SaveResult(out1$vs,"out1.txt")
out2=simulate(L,n,beta,"B")
SaveResult(out2$vs,"out2.txt")
out3=simulate(L,n,beta,"C")
SaveResult(out3$vs,"out3.txt")
out4=simulate(L,n,beta,"D",penalty1="log")
SaveResult(out4$vs,"out4.txt")

out4_500=simulate(L,n,beta,"D")
SaveResult(out4_500$vs,"out4_500.txt")
out3_500=simulate(L,n,beta,"C")
SaveResult(out3_500$vs,"out3_500.txt")
out1_500=simulate(L,n,beta,"A")
SaveResult(out1_500$vs,"out1_500.txt")
out2_500=simulate(L,n,beta,"B")
SaveResult(out2_500$vs,"out2_500.txt")




outADL1_500=simulate(L,n,beta,"A",method="ADL")
SaveResult(outADL1_500$vs,"outADL1_500.txt")
outADL2_500=simulate(L,n,beta,"B",method="ADL")
SaveResult(outADL2_500$vs,"outADL2_500.txt")
outADL3_500=simulate(L,n,beta,"C",method="ADL")
SaveResult(outADL3_500$vs,"outADL3_500.txt")
outADL4_500=simulate(L,n,beta,"D",method="ADL")
SaveResult(outADL4_500$vs,"outADL4_500.txt")


out1_500=out1
out2_500=out2
out3_500=out3
out4_500=out4
#-------------------t test for symmetricize--------------------#
t=20
L=50
n=50
p=8
beta1=c(3,2,1.5,1,1,1,1,1)
beta2=c(3,2,1.5,0,0,0,0,0)
beta=beta2
rate1=rep(0,t)
rate2=rep(0,t)
for(i in 1:t)
{
  out3=simulate(L,n,beta,"C",seed=i)
  rate2[i]=out3$vs[1]
  #out3=simulate(L,n,beta,"C3",seed=i)
  #rate2[i]=out3$vs[1]
}
t.test(rate2-rate1)

#-------------real data analysis-------------#
require(MASS)
y=log(Boston$medv)
x=Boston[,-14]
out=list(y=y,x=x)

out=PrepareData(screenNum=1000)
p=dim(out$x)[2]
n=dim(out$x)[1]
init=sparseLTS(out$x,out$y)
beta0=SetBeta0(init$coefficients)
w0=ifelse(init$wt==1,0.99,0.01)

beta0=rep(1,p+1)
w0=rep(0.99,n)
nlambda1=50
nlambda2=100
res=srcdreg(out$x,out$y,nlambda1=50,nlambda2=100,beta0=beta0,w0=w0,delta=0.000001,maxIter=1000)

res=rcdreg(out$x,out$y,penalty="ADL",nlambda1=nlambda1,nlambda2=nlambda2,beta0=beta0,w0=w0,delta=0.000001,maxIter=1000)



#---------ADL----------#
init2=InitParam(out$x,out$y,method="LASSO")''
beta0=SetBeta0(init$beta)
w0=ifelse(init$weight==1,(0.99+max(init$weight[init$weight!=1])*0.01),init$weight)
res=rcdreg(out$x,out$y,penalty="ADL",nlambda1=2,nlambda2=nlambda2,beta0=beta0,w0=w0,delta=0.000001,maxIter=1000)
index=BIC(res$wloss[1,],apply(matrix(res$beta[1,,],100,p)!=0+0,1,sum),n,p)
bbb=res$beta[1,index,]

ltsReg(out$x,out$y)


res=slim(out$x,out$y,method="lq",q=1,nlambda=100,verbose=FALSE)
res=InitParam(out$x,out$y,method="LAD")



#============Boston Housing==========#
require(MASS)
require(robustHD)
y=log(Boston$medv)
x=Boston[,-14]
out=list(y=y,x=x)
p=dim(out$x)[2]
n=dim(out$x)[1]
init=sparseLTS(out$x,out$y)
beta0=SetBeta0(init$coefficients)
w0=ifelse(init$wt==1,0.99,0.01)
nlambda1=50
nlambda2=100
# first pwls-vs
res=srcdreg(out$x,out$y,nlambda1=nlambda1,nlambda2=nlambda2,beta0=beta0,w0=w0,delta=0.000001,maxIter=1000,updateInitial=FALSE)

#pwls under selected space of beta
select=(res$beta!=0)
xx=x[,select]
res2=pwlsreg(xx,out$y,nlambda=50,w0=res$w,delta=0.0001,maxIter=1000)



# calculate varinace of beta hat

#=========end of Boston Housing=====#

#-----------------------------------------------------------------------------------------------------------------

#======real data analysis============#



#======end of data analysis===========#