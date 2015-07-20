#n=50,p=8
L=100
n=50
p=8
beta=c(3,2,1.5,0,0,0,0,0)

#MMNNG collect data- low dimension
simulate(L,n,beta,"A",method="MMNNG_DATA",seed=NULL)
simulate(L,n,beta,"B",method="MMNNG_DATA",seed=NULL)
simulate(L,n,beta,"C",method="MMNNG_DATA",seed=NULL)
simulate(L,n,beta,"D",method="MMNNG_DATA",seed=NULL)
simulate(L,n,beta,"D2",method="MMNNG_DATA",seed=NULL)


#MMNNG-low dimension
outA_MM50=simulate(L,n,beta,"A",method="MMNNG",useDataFile=TRUE)
outB_MM50=simulate(L,n,beta,"B",method="MMNNG",useDataFile=TRUE)
outC_MM50=simulate(L,n,beta,"C",method="MMNNG",useDataFile=TRUE)
outD_MM50=simulate(L,n,beta,"D",method="MMNNG",useDataFile=TRUE)
outD2_MM50=simulate(L,n,beta,"D2",method="MMNNG",useDataFile=TRUE)

#LTS-low dimension
outA_LTS50=simulate(L,n,beta,"A",method="LTS",useDataFile=TRUE)
outB_LTS50=simulate(L,n,beta,"B",method="LTS",useDataFile=TRUE)
outC_LTS50=simulate(L,n,beta,"C",method="LTS",useDataFile=TRUE)
outD_LTS50=simulate(L,n,beta,"D",method="LTS",useDataFile=TRUE)
outD2_LTS50=simulate(L,n,beta,"D2",method="LTS",useDataFile=TRUE)

#PWLQ-low dimension
confirm=simulate(2,n,beta,"A",method="PWLQ",initial="LTS",seed=NULL,useDataFile=TRUE)
outA_50=simulate(L,n,beta,"A",method="PWLQ",initial="LTS",seed=NULL,useDataFile=TRUE)
outB_50=simulate(L,n,beta,"B",method="PWLQ",initial="LTS",seed=NULL,useDataFile=TRUE)
outC_50=simulate(L,n,beta,"C",method="PWLQ",initial="LTS",seed=NULL,useDataFile=TRUE)
outD2_50=simulate(L,n,beta,"D2",method="PWLQ",initial="LTS",seed=NULL,useDataFile=TRUE)
#ROSS
matLabDir="D:\\matlab\\Sim"
matlab=PrepareMatlab(matLabDir)
outA_ROSS50=simulate(L,n,beta,"A",method="ROSS",matlab=matlab,useDataFile=TRUE)
outB_ROSS50=simulate(L,n,beta,"B",method="ROSS",matlab=matlab,useDataFile=TRUE)
outC_ROSS50=simulate(L,n,beta,"C",method="ROSS",matlab=matlab,useDataFile=TRUE)
outD_ROSS50=simulate(L,n,beta,"D",method="ROSS",matlab=matlab,useDataFile=TRUE)
outD2_ROSS50=simulate(L,n,beta,"D2",method="ROSS",matlab=matlab,useDataFile=TRUE)
#ADL
outA_ADL50=simulate(L,n,beta,"A",method="ADL",useDataFile=TRUE)
outB_ADL50=simulate(L,n,beta,"B",method="ADL",useDataFile=TRUE)
outC_ADL50=simulate(L,n,beta,"C",method="ADL",useDataFile=TRUE)
outD_ADL50=simulate(L,n,beta,"D",method="ADL",useDataFile=TRUE)
outD2_ADL50=simulate(L,n,beta,"D2",method="ADL",useDataFile=TRUE)
#-------------------------------------------------------------------------------------

#n=100,500
L=100
n=100
p=500
beta=c(rep(2,10),rep(0,p-10))

#MMNNG collect data- high dimension
simulate(L,n,beta,"A",method="MMNNG_DATA",seed=NULL)
simulate(L,n,beta,"B",method="MMNNG_DATA",seed=NULL)
simulate(L,n,beta,"C",method="MMNNG_DATA",seed=NULL)
simulate(L,n,beta,"D",method="MMNNG_DATA",seed=NULL)

#ROSS
matLabDir="D:\\matlab\\Sim"
matlab=PrepareMatlab(matLabDir)
outA_ROSS500=simulate(L,n,beta,model="D",method="ROSS",matlab=matlab)
outA_ROSS500=simulate(L,n,beta,model="C",method="ROSS",matlab=matlab)
outA_ROSS500=simulate(L,n,beta,model="B",method="ROSS",matlab=matlab)
outA_ROSS500=simulate(L,n,beta,model="A",method="ROSS",matlab=matlab)
#pwlq
outA_500=simulate(L,n,beta,"A",method="PWLQ",initial="plain",seed=2015,updateInitial=TRUE)
outC_500=simulate(L,n,beta,"C",method="PWLQ",initial="plain",seed=2015,updateInitial=TRUE)
outD_500=simulate(L,n,beta,"D",method="PWLQ",initial="plain",seed=2015,updateInitial=TRUE)
outB_500=simulate(L,n,beta,"B",method="PWLQ",initial="plain",seed=2015,updateInitial=TRUE)
outD2_500=simulate(L,n,beta,"D2",method="PWLQ",initial="plain",seed=2015,updateInitial=TRUE)
outC_500=simulate(L,n,beta,"C",method="PWLQ",initial="plain",seed=2015,updateInitial=TRUE)
#ADL
outA_ADL500=simulate(L,n,beta,"A",method="ADL",seed=2015)
outC_ADL500=simulate(L,n,beta,"C",method="ADL",seed=2015)
outB_ADL500=simulate(L,n,beta,"B",method="ADL",useDataFile=TRUE)
outD_ADL500=simulate(L,n,beta,"D",method="ADL",useDataFile=TRUE) 
outD2_ADL500=simulate(L,n,beta,"D",method="ADL",seed=2015)
#LTS
SaveResult=function(res,file,dir="C:\\Users\\Administrator\\Dropbox\\result\\")
{
  write.table(res,paste(dir,file,sep=""))
}
outD2_LTS500=simulate(L,n,beta,"D2",method="LTS",seed=2015)
outD2_500_test2=simulate(L,n,beta,"D2",method="PWLQ",initial="LTS",seed=2015,updateInitial=TRUE)
SaveResult(outD2_500$vs,"outD2_500.txt")
outA_500_test=simulate(L,n,beta,"A",method="PWLQ",initial="plain",seed=2015,updateInitial=TRUE,updateInitialTimes=100)
SaveResult(outA_500_test$vs,"outA_500_test.txt")
#-------------------------------example 2---------------------#

#high n=100,p=1000
n=100
p=1000
beta=rep(0,p)
beta[c(1,7)]=1.5
beta[2]=0.5
beta[c(4,11)]=1
res=simulate(1,n,beta,"HA",method="PWLQ",initial="LTS",seed=2015)

#low
n=150
p=50
k=6
L=100
beta=c(rep(1,k),rep(0,p-k))
#MMNNG collect data
simulate(L,n,beta,"LA",method="MMNNG_DATA",seed=NULL)
simulate(L,n,beta,"LB",method="MMNNG_DATA",seed=NULL)
simulate(L,n,beta,"LC",method="MMNNG_DATA",seed=NULL)
simulate(L,n,beta,"LD",method="MMNNG_DATA",seed=NULL)

#PWLQ low dimension
outLA_PWLQ=simulate(L,n,beta=beta,model="LA",method="PWLQ",initial="LTS",useDataFile=TRUE)
outLB_PWLQ=simulate(L,n,beta=beta,model="LB",method="PWLQ",initial="LTS",useDataFile=TRUE)
outLc_PWLQ=simulate(L,n,beta=beta,model="LC",method="PWLQ",initial="LTS",seed=2015)

#Spare LTS low dimension
outLA_LTS=simulate(L,n,beta=beta,model="LA",method="LTS",useDataFile=TRUE)
outLB_LTS=simulate(L,n,beta=beta,model="LB",method="LTS",useDataFile=TRUE)
outLc_LTS=simulate(L,n,beta=beta,model="LC",method="LTS",seed=2015)

#MMNNG-low dimension
outLA_MM=simulate(L,n,beta,"LA",method="MMNNG",useDataFile=TRUE)
outLB_MM=simulate(L,n,beta,"LB",method="MMNNG",useDataFile=TRUE)
#outC_MM50=simulate(L,n,beta,"C",method="MMNNG",useDataFile=TRUE)
#outD_MM50=simulate(L,n,beta,"D",method="MMNNG",useDataFile=TRUE)

#ADL low dimension
outLA_ADL=simulate(L,n,beta,"LA",method="ADL",useDataFile=TRUE)
outLB_ADL=simulate(L,n,beta,"LB",method="ADL",useDataFile=TRUE)

out=GenerateDataSep(n=10,p=10,k=2)
#----load data------#
load("outB_500.rda")
load("outB_ADL500.rda")
load("outB_LTS500.rda")

#======real data analysis(air pollution)============#
airRaw=read.delim("airPollution.txt")
airData=airRaw[-21,-1]
airData$HCPot=log(airData$HCPot)
airData$NOxPot=log(airData$NOxPot)
airData$S02Pot=log(airData$S02Pot)
tempData=t(airData)-apply(airData,2,median)
data=t(tempData/apply(abs(tempData),1,median))
y=data[,5]
x=data[,-c(5)]
out=list(y=y,x=x)
colnames=names(airData[,-c(5)])

#OLS
lm1=lm(y~x)
lm0=lm(as.vector(airData[,5])~as.matrix(airData[,-c(5)]))
# pwls-vs
p=dim(out$x)[2]
n=dim(out$x)[1]
require(robustHD)
init=sparseLTS(out$x,out$y)
beta0=SetBeta0(init$coefficients)
w0=ifelse(init$wt==1,0.99,0.01)
nlambda1=50
nlambda2=100
res=srcdreg(out$x,out$y,penalty1="1-w0",nlambda1=nlambda1,nlambda2=nlambda2,beta0=beta0,w0=w0,delta=0.000001,maxIter=1000,
            intercept=TRUE,standardize=FALSE,updateInitial=FALSE,criterion="BIC")
colnames[res$beta[-1]!=0]
studres=studres(lm1)
plot(studres)
abline(2.5,0)
identify(studres)
#MMNNG
source("mmnngreg.R")
res_MM=mmnngreg(as.matrix(out$x),out$y)
res_MM$betac[-1]=colnames

#======real data analysis(NCI-60)============#
#read data
nci_pro=read.delim("nci60_Protein__Lysate_Array_log2.txt")
nci_gene0=read.delim("RNA__Affy_HG_U133(A_B)_GCRMA.txt")
nci_temp=read.delim("GPL96-15653.txt")
nci_gene=nci_gene0[nci_gene0$Probe.id..b.%in%nci_temp$ID,]
#obtain X and y

x=t(nci_gene[,-c(1:9,49,70)])
x=matrix(x,nrow=dim(x)[1])
y=as.vector(nci_pro[92,-c(1:4,44,65)])
y=as.vector(matrix(data=y,ncol=1))
out=list(y=y,x=x)

#
require(robustHD)
class(out$x)<-"numeric"
res=sparseLTS(out$x,out$y)
#compute
# pwls-vs
p=dim(out$x)[2]
n=dim(out$x)[1]
beta0=rep(1,p+1)
w0=rep(0.99,n)
nlambda1=50
nlambda2=100
res=srcdreg(out$x,out$y,penalty1="1-w0",nlambda1=nlambda1,nlambda2=nlambda2,beta0=beta0,w0=w0,delta=0.000001,maxIter=1000,
            intercept=TRUE,standardize=FALSE,updateInitial=FALSE,criterion="BIC")
colnames[res$beta[-1]!=0]
#======end of data analysis===========#