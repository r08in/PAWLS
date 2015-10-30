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

#PAWLS-low dimension
confirm=simulate(L,n,beta,"A",method="PAWLS",initial="LTS",seed=NULL,useDataFile=TRUE)
outA_50=simulate(L,n,beta,"A",method="PAWLS",initial="LTS",seed=NULL,useDataFile=TRUE)
outB_50=simulate(L,n,beta,"B",method="PAWLS",initial="LTS",seed=NULL,useDataFile=TRUE)
outC_50=simulate(L,n,beta,"C",method="PAWLS",initial="LTS",seed=NULL,useDataFile=TRUE)
outD_50=simulate(L,n,beta,"D",method="PAWLS",initial="LTS",seed=NULL,useDataFile=TRUE)
outD2_50=simulate(L,n,beta,"D2",method="PAWLS",initial="LTS",seed=NULL,useDataFile=TRUE)
#ROSS
matLabDir="D:\\matlab\\Sim"
matlab=PrepareMatlab(matLabDir)
outA_ROSS50=simulate(L,n,beta,"A",method="ROSS",matlab=matlab,useDataFile=TRUE)
outB_ROSS50=simulate(L,n,beta,"B",method="ROSS",matlab=matlab,useDataFile=TRUE)
outC_ROSS50=simulate(L,n,beta,"C",method="ROSS",matlab=matlab,useDataFile=TRUE)
outD_ROSS50=simulate(L,n,beta,"D",method="ROSS",matlab=matlab,useDataFile=TRUE)
outD2_ROSS50=simulate(L,n,beta,"D2",method="ROSS",matlab=matlab,useDataFile=TRUE)
#ADL
confirm_ADL50=simulate(L,n,beta,"A",method="ADL",useDataFile=TRUE)
outA_ADL50_adl=simulate(L,n,beta,"A",method="ADL",useDataFile=TRUE)
outB_ADL50=simulate(L,n,beta,"B",method="ADL",useDataFile=TRUE)
outC_ADL50=simulate(L,n,beta,"C",method="ADL",useDataFile=TRUE)
outD_ADL50=simulate(L,n,beta,"D",method="ADL",useDataFile=TRUE)
outD2_ADL50=simulate(L,n,beta,"D2",method="ADL",useDataFile=TRUE)
#-------------------------------------------------------------------------------------

#n=100,500
L=100
n=100
p=500
num=10
beta=c(rep(2,num),rep(0,p-num))

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
#PAWLS
confrim_500=simulate(L,n,beta,"A",method="PAWLS",seed=2015,updateInitial=TRUE)
outA_500=simulate(L,n,beta,"A",method="PAWLS",initial="uniform",seed=2015,updateInitial=TRUE)
outC_500=simulate(L,n,beta,"C",method="PAWLS",initial="uniform",seed=2015,updateInitial=TRUE)
outD_500=simulate(L,n,beta,"D",method="PAWLS",initial="uniform",seed=2015,updateInitial=TRUE)
outB_500=simulate(L,n,beta,"B",method="PAWLS",initial="uniform",seed=2015,updateInitial=TRUE)
outD2_500=simulate(L,n,beta,"D2",method="PAWLS",initial="uniform",seed=2015,updateInitial=TRUE)
outD3_500=simulate(L,n,beta,"D3",method="PAWLS",initial="uniform",seed=2015,updateInitial=TRUE)
outC_500=simulate(L,n,beta,"C",method="PAWLS",initial="uniform",seed=2015,updateInitial=TRUE)
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
outA_LTS500=simulate(L,n,beta,"A",method="LTS",seed=2015)
outB_LTS500=simulate(L,n,beta,"B",method="LTS",seed=2015)
outC_LTS500=simulate(L,n,beta,"C",method="LTS",seed=2015)
outD2_LTS500=simulate(L,n,beta,"D2",method="LTS",seed=2015)


outD2_500_test2=simulate(L,n,beta,"D2",method="PAWLS",initial="LTS",seed=2015,updateInitial=TRUE)
SaveResult(outD2_500$vs,"outD2_500.txt")
outA_500_test=simulate(L,n,beta,"A",method="PAWLS",initial="uniform",seed=2015,updateInitial=TRUE,updateInitialTimes=100)
SaveResult(outA_500_test$vs,"outA_500_test.txt")
#-------------------------------example 2---------------------#

#high n=100,p=1000
n=100
p=1000
beta=rep(0,p)
beta[c(1,7)]=1.5
beta[2]=0.5
beta[c(4,11)]=1
res=simulate(1,n,beta,"HA",method="PAWLS",initial="LTS",seed=2015)

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

#PAWLS low dimension
outLA_PAWLS=simulate(L,n,beta=beta,model="LA",method="PAWLS",initial="LTS",useDataFile=TRUE)
outLB_PAWLS=simulate(L,n,beta=beta,model="LB",method="PAWLS",initial="LTS",useDataFile=TRUE)
outLc_PAWLS=simulate(L,n,beta=beta,model="LC",method="PAWLS",initial="LTS",seed=2015)

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
x=data[,-c(5,16)]
out=list(y=y,x=x)
colnames=names(airData[,-c(5,16)])
colnames[7]="NonWhite"
#OLS
lm1=lm(y~x)
lm0=lm(as.vector(airData[,5])~as.matrix(airData[,-c(5)]))

# pwls-vs
res_air=srcdreg(out$x,out$y,initial="LTS",search="crossDynamic")
res2_air=srcdreg(out$x,out$y,initial="LTS",search="all")
colnames[res_air$beta[-1]!=0]
studres=studres(lm1)
plot(studres)
abline(2.5,0)
identify(studres)
paramPlot(res2_air,xlimw=c(6,12))
#getlabel
label_beta=rep(colnames,100)
paramPlot(res2_air,label2=label_beta)

#ADL
res_air_ADL=srcdreg(out$x,out$y,nlambda1=2,initial="LAD",search="fixw")
colnames[res_air_ADL$beta[-1]!=0]

#sLTS
require(robustHD)
res_air_LTS=sparseLTS(out$x,out$y,intercept=TRUE)
colnames[res_air_LTS$coefficients[-1]!=0]

# MMNNG
source("mmnngreg.R")
res_air_MMNNG=mmnngreg(out$x,out$y)
colnames[res_air_MMNNG$betac[-1]!=0]

#SROS
x=AddIntercept(out$x)
setVariable(matlab, X=x)
setVariable(matlab, y=out$y)
evaluate(matlab,"[betaRoss]=RossSimulate(X,y)")
res_air_Ross=getVariable(matlab, "betaRoss")
colnames[res_air_Ross$betaRoss[-1]!=0]

#======real data analysis(NCI-60)============#
#read data
nci_pro=read.delim("nci60_Protein__Lysate_Array_log2.txt")
nci_gene0=read.delim("RNA__Affy_HG_U133(A_B)_GCRMA.txt")
nci_temp=read.delim("GPL96-15653.txt")
nci_gene=nci_gene0[nci_gene0$Probe.id..b.%in%nci_temp$ID,]
cnames=nci_gene[,2]
#obtain X and y

x=t(nci_gene[,-c(1:9,49,70)])
x=matrix(as.numeric(x),nrow=dim(x)[1])
colnames(x)<-cnames
y=as.numeric(nci_pro[92,-c(1:4,44,65)])
screenNum=500
out=ScreenData(y,x,screenNum)
colname=colnames(out$x)

#ADL
res_nci_ADL=srcdreg(out$x,out$y,nlambda1=2,initial="LAD",search="fixw")
colnames[res_air_ADL$beta[-1]!=0]

#LTS
require(robustHD)
#class(out$x)<-"numeric"
res=sparseLTS(out$x,out$y)
res$coefficients[res$coefficients!=0]

#a-lasso
res_adl_nci=srcdreg(out$x,out$y,nlambda1=2,initial="LAD",search="fixw")
names(res_adl_nci$beta)<-c("intercept",colname)
res_adl_nci$beta[res_adl_nci$beta!=0]

# pwls-vs
res_nci=srcdreg(out$x,out$y,nlambda2=100,updateInitialTimes=4,search="crossDynamic",criterion="BIC") #crossDynamic
res_nci_AIC=res_nci
res_nci_BIC=res_nci
names(res_nci_AIC$beta)<-c("intercept",colname)
names(res_nci_BIC$beta)<-c("intercept",colname)
res_nci_AIC$beta[res_nci_AIC$beta!=0]
res_nci_BIC$beta[res_nci_BIC$beta!=0]
res2_nci=srcdreg(out$x,out$y,beta0=res_nci$beta0,w0=res_nci$w0,
                 lambda1=res_nci$lambda1s,lambda2=res_nci$lambda2s,
                  updateInitialTimes=4,search="all",criterion="AIC")
paramPlot(res2_nci,param="both",xlimb=c(0,11),xlimw=c(5.5,11))
paramPlot(res2_nci,param="both",xlimw=c(5.5,11))
names(res2_nci$beta)<-c("intercept",colname)
res2_nci$beta[res2_nci$beta!=0]

colnames[res$beta[-1]!=0]

#seq(0.02013544,0,length.out=100)
corr=cor(out$x[,1],out$x[,-1])
View(corr[order(abs(corr),decreasing=TRUE)])

res_nci2=srcdreg(out$x[,-1],out$y-(out$x[,1]*res_nci$beta[2]),updateInitialTimes=4) #1

i=index[res2_nci$beta[-1]!=0]
res2_nci$beta[res2_nci$beta!=0]
colname[i]
#======end of data analysis===========#

#=======temp==========#
View(outD2_50$beta)
View(outD2_LTS50$beta)
se=c(outD2_50$se,outD2_LTS50$se)
group=c(rep("PAWLS",100),rep("SparseLTS",100))
m=data.frame(se=se,group=group)
boxplot(se~group,data=m,main="model D",outline=FALSE)


#=======end temp=======#