#----------test between 1-wo and log wo---------------#
#case when 1-wo
t=20
L=50
n=50
p=8
beta1=c(3,2,1.5,1,1,1,1,1)
beta2=c(3,2,1.5,0,0,0,0,0)
beta=beta2
rate1=matrix(nrow=4,ncol=t,0) #1-w0
rate2=matrix(nrow=4,ncol=t,0) #log
rate3=matrix(nrow=4,ncol=t,0) #standardize=TRUE
rate4=matrix(nrow=4,ncol=t,0) #updateInitial=FALSE
for(i in 1:t) # for 1-wo
{
  rate1[1,i]=simulate(L,n,beta,"A",seed=i)$vs[1]
  rate1[2,i]=simulate(L,n,beta,"B",seed=i)$vs[1]
  rate1[3,i]=simulate(L,n,beta,"C",seed=i)$vs[1]
  rate1[4,i]=simulate(L,n,beta,"D",seed=i)$vs[1]
}
for(i in 1:t)#for log wo
{
  rate2[1,i]=simulate(L,n,beta,"A",seed=i)$vs[1]
  rate2[2,i]=simulate(L,n,beta,"B",seed=i)$vs[1]
  rate2[3,i]=simulate(L,n,beta,"C",seed=i)$vs[1]
  rate2[4,i]=simulate(L,n,beta,"D",seed=i)$vs[1]
}
for(i in 1:t)#for standardize=TRUE
{
  rate3[1,i]=simulate(L,n,beta,"A",seed=i,penalty1="1-w0",standardize=TRUE)$vs[1]
  rate3[2,i]=simulate(L,n,beta,"B",seed=i,penalty1="1-w0",standardize=TRUE)$vs[1]
  rate3[3,i]=simulate(L,n,beta,"C",seed=i,penalty1="1-w0",standardize=TRUE)$vs[1]
  rate3[4,i]=simulate(L,n,beta,"D",seed=i,penalty1="1-w0",standardize=TRUE)$vs[1]
}
for(i in 1:t)#for updateInitial=FALSE
{
  rate4[1,i]=simulate(L,n,beta,"A",seed=i,updateInitial=FALSE)$vs[1]
  rate4[2,i]=simulate(L,n,beta,"B",seed=i,updateInitial=FALSE)$vs[1]
  rate4[3,i]=simulate(L,n,beta,"C",seed=i,updateInitial=FALSE)$vs[1]
  rate4[4,i]=simulate(L,n,beta,"D",seed=i,updateInitial=FALSE)$vs[1]
}
apply(rate1-rate2,1,t.test)


#----------end test------------------------------------#wkzz    

#----------find possible seed for MM-nng in low dimension--------#
L=100
n=50
p=8
beta1=c(3,2,1.5,1,1,1,1,1)
beta2=c(3,2,1.5,0,0,0,0,0)
beta=beta2
beta0=matrix(nrow=L,ncol=p,0)
w0=matrix(nrow=L,ncol=n,0)
beta_MM=matrix(nrow=L,ncol=p,0)

t=which(apply(beta0,1,sum)==0)[1]
for(i in t:L)
{
  res=simulate(1,n,beta,"D",seed=NULL,updateInitial=FALSE,initial="LTS",method="COMPARE")
  write.table(res$beta,"beta.txt",append=TRUE)
  write.table(res$w,"w.txt",append=TRUE)
  write.table(res$beta1,"beta_MMNNG.txt",append=TRUE)
}

#----------end finding--------------------------------------------#

#----------data collection-----------#


#n=50,p=8
L=100
n=50
p=8
beta1=c(3,2,1.5,1,1,1,1,1)
beta2=c(3,2,1.5,0,0,0,0,0)
beta=beta2

#MMNNG collect data- low dimension
simulate(L,n,beta,"A",method="MMNNG_DATA",seed=NULL)
simulate(L,n,beta,"B",method="MMNNG_DATA",seed=NULL)
simulate(L,n,beta,"C",method="MMNNG_DATA",seed=NULL)
simulate(L,n,beta,"D",method="MMNNG_DATA",seed=NULL)


#MMNNG-low dimension
outA_MM50=simulate(L,n,beta,"A",method="MMNNG",useDataFile=TRUE)
outB_MM50=simulate(L,n,beta,"B",method="MMNNG",useDataFile=TRUE)
outC_MM50=simulate(L,n,beta,"C",method="MMNNG",useDataFile=TRUE)
outD_MM50=simulate(L,n,beta,"D",method="MMNNG",useDataFile=TRUE)
#LTS- dimension
outA_LTS50=simulate(L,n,beta,"A",method="LTS",useDataFile=TRUE)
outB_LTS50=simulate(L,n,beta,"B",method="LTS",useDataFile=TRUE)
outC_LTS50=simulate(L,n,beta,"C",method="LTS",useDataFile=TRUE)
outD_LTS50=simulate(L,n,beta,"D",method="LTS",useDataFile=TRUE)
#PWLQ-low dimension
outA_50=simulate(L,n,beta,"A",method="PWLQ",initial="LTS",seed=NULL,useDataFile=TRUE)
outB_50=simulate(L,n,beta,"B",method="PWLQ",initial="LTS",seed=NULL,useDataFile=TRUE)
outC_50=simulate(L,n,beta,"C",method="PWLQ",initial="LTS",seed=NULL,useDataFile=TRUE)
outD_50=simulate(L,n,beta,"D",method="PWLQ",initial="LTS",seed=NULL,useDataFile=TRUE)

#ROSS
matLabDir="D:\\matlab\\Sim"
matlab=PrepareMatlab(matLabDir)
outA_ROSS50=simulate(L,n,beta,"A",method="ROSS",matlab=matlab,useDataFile=TRUE)
outB_ROSS50=simulate(L,n,beta,"B",method="ROSS",matlab=matlab,useDataFile=TRUE)
outC_ROSS50=simulate(L,n,beta,"C",method="ROSS",matlab=matlab,useDataFile=TRUE)
outD_ROSS50=simulate(L,n,beta,"D",method="ROSS",matlab=matlab,useDataFile=TRUE)
#ADL
outA_ADL50=simulate(L,n,beta,"A",method="ADL",useDataFile=TRUE)
outB_ADL50=simulate(L,n,beta,"B",method="ADL",useDataFile=TRUE)
outC_ADL50=simulate(L,n,beta,"C",method="ADL",useDataFile=TRUE)
outD_ADL50=simulate(L,n,beta,"D",method="ADL",useDataFile=TRUE)
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
outA_500_test=simulate(L,n,beta,"A",method="PWLQ",initial="plain",seed=2015,updateInitial=TRUE)
outC_500=simulate(L,n,beta,"C",method="PWLQ",initial="plain",seed=2015,updateInitial=TRUE)
outD_500=simulate(L,n,beta,"D",method="PWLQ",initial="plain",seed=2015,updateInitial=TRUE)
outB_500=simulate(L,n,beta,"B",method="PWLQ",initial="plain",seed=2015,updateInitial=TRUE)

#ADL
outA_ADL500=simulate(L,n,beta,"A",method="ADL",seed=2015)
outC_ADL500=simulate(L,n,beta,"C",method="ADL",seed=2015)
outB_ADL500=simulate(L,n,beta,"B",method="ADL",useDataFile=TRUE)
outD_ADL500=simulate(L,n,beta,"D",method="ADL",useDataFile=TRUE)

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
