#RRMM
L=100
n=50
p=8
#beta=rep(1,p)
beta=c(3,2,1,1,1,1,0.1,0.1)
matLabDir="D:\\matlab\\RRMM"
matlab=PrepareMatlab(matLabDir)


#RRMM
rA_RRMM50=simulate2(L,n,beta,"RA",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2016)
rB_RRMM50=simulate2(L,n,beta,"RB",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2016)
rC_RRMM50=simulate2(L,n,beta,"RC",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2016)
rD_RRMM50=simulate2(L,n,beta,"RD",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2016)

#RRREG
outrA_50=simulate2(L,n,beta,"RA",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2016,
                   lambda1=seq(from=0.001,to=1,length.out=50),lambda2=seq(from=0.001,to=0.1,length.out=30))
outrB_50=simulate2(L,n,beta,"RB",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",,seed=2016,
                   lambda1=seq(from=0.11,to=5,length.out=100),lambda2=seq(from=0.005,to=0.5,length.out=100))
outrc_50=simulate2(L,n,beta,"RC",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",,seed=2016,
                   lambda1=seq(from=1,to=2,length.out=100),lambda2=seq(from=0.005,to=0.5,length.out=100))
outrD_50=simulate2(L,n,beta,"RD",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",,seed=2016,
                   lambda1=seq(from=0.001,to=1,length.out=50),lambda2=seq(from=0.005,to=0.5,length.out=100))

#A-RRREG

#test
tA_50=simulate2(L,n,beta,"RA",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2016,
                   lambda1=seq(from=0.2,to=2,length.out=50),lambda2=seq(from=0,to=0,length.out=2))
tB_50=simulate2(L,n,beta,"RB",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",,seed=2016,
                   lambda1=seq(from=0.11,to=5,length.out=100),lambda2=seq(from=0,to=0,length.out=2))
tc_50=simulate2(L,n,beta,"RC",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",,seed=2016,
                   lambda1=seq(from=1,to=2,length.out=100),lambda2=seq(from=0,to=0,length.out=2))
tD_50=simulate2(L,n,beta,"RD",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",,seed=2016,
                   lambda1=seq(from=0.001,to=1,length.out=50),lambda2=seq(from=0,to=0,length.out=2))
set.seed(2016)
out=GenerateDataByModel(n,beta=beta,model="A",dataType="Ridge")
res=rrreg(x=out$x,y=out$y,lambda1=seq(from=0.001,to=1,length.out=30),lambda2=seq(from=0.01,to=2,length.out=30),intercept=TRUE)
GetRobustPe(AddIntercept(out$x),y=out$y,betaHat = res$beta)
res=rrreg(x=out$x,y=out$y,lambda1=seq(from=10,to=20,length.out=100),lambda2=seq(from=0.0001,to=1,length.out=100),intercept="false")

outC=GenerateDataByModel(n,beta=beta,model="A",r=0,dataType="Ridge")

res=rrreg(x=outC$x,y=outC$y,lambda1=1,lambda2=0,intercept="false")


