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