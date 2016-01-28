#RRMM
L=10
n=50
p=8
beta=rep(1,p)
#beta=c(3,2,1.5,1.5,1,1,0.5,0.5)
matLabDir="D:\\matlab\\RRMM"
matlab=PrepareMatlab(matLabDir)


#RRMM
outA_RRMM50=simulate2(L,n,beta,"A",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2016)
outB_RRMM50=simulate2(L,n,beta,"B",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2016)
outC_RRMM50=simulate2(L,n,beta,"C",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2016)
outD2_RRMM50=simulate2(L,n,beta,"D2",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2016)

#RRREG
outrA_50=simulate2(L,n,beta,"A",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2016,
                   lambda1=seq(from=0.001,to=1,length.out=50),lambda2=seq(from=0.001,to=0.1,length.out=30))
outrc_50=simulate2(L,n,beta,"C",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",,seed=2016,
                   lambda1=seq(from=1,to=2,length.out=100),lambda2=seq(from=0.005,to=0.5,length.out=100))
outrB_50=simulate2(L,n,beta,"B",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",,seed=2016,
                   lambda1=seq(from=0.11,to=5,length.out=100),lambda2=seq(from=0.005,to=0.5,length.out=100))
outrD_50=simulate2(L,n,beta,"D",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",,seed=2016,
                   lambda1=seq(from=0.001,to=1,length.out=50),lambda2=seq(from=0.005,to=0.5,length.out=100))

#A-RRREG


set.seed(2016)
out=GenerateDataByModel(n,beta=beta,model="A",dataType="Ridge")
res=rrreg(x=out$x,y=out$y,lambda1=seq(from=0.001,to=1,length.out=30),lambda2=seq(from=0.01,to=2,length.out=30),intercept=TRUE)
GetRobustPe(AddIntercept(out$x),y=out$y,betaHat = res$beta)
res=rrreg(x=out$x,y=out$y,lambda1=seq(from=10,to=20,length.out=100),lambda2=seq(from=0.0001,to=1,length.out=100),intercept="false")

outC=GenerateDataByModel(n,beta=beta,model="A",r=0,dataType="Ridge")

res=rrreg(x=outC$x,y=outC$y,lambda1=1,lambda2=0,intercept="false")


