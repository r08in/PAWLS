#----------------------------low dimension-----------------------------------#
#RRMM
L=100
n=50
p=8
#beta=rep(1,p)
beta=c(3,2,1,1,1,1,0,0)
#beta=rep(0,p)

matLabDir="D:\\matlab\\RRMM"
matlab=PrepareMatlab(matLabDir)

#RRMM50
rA_RRMM50=simulate2(L,n,beta,"RA",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017)
rB_RRMM50=simulate2(L,n,beta,"RB",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017)
rC_RRMM50=simulate2(L,n,beta,"RC",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017)
rD_RRMM50=simulate2(L,n,beta,"RD",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017)

#RRMM50.3
rA_RRMM50.3=simulate2(L,n,beta,"RA",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017,pro=0.3)
rB_RRMM50.3_test=simulate2(L,n,beta,"RB",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017,pro=0.3)
rC_RRMM50.3=simulate2(L,n,beta,"RC",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017,pro=0.3)
rD_RRMM50.3=simulate2(L,n,beta,"RD",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017,pro=0.3)

#RRMM1000
rA_RRMM1000=simulate2(L,n,beta,"RA",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017)
rB_RRMM1000=simulate2(L,n,beta,"RB",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017)
rC_RRMM1000=simulate2(L,n,beta,"RC",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017)
rD_RRMM1000=simulate2(L,n,beta,"RD",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017)


#RRMM1000.3
rA_RRMM1000.3=simulate2(L,n,beta,"RA",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017,pro=0.3)
rB_RRMM1000.3=simulate2(L,n,beta,"RB",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017,pro=0.3)
rC_RRMM1000.3=simulate2(L,n,beta,"RC",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017,pro=0.3)
rD_RRMM1000.3=simulate2(L,n,beta,"RD",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017,pro=0.3)

#RRREG50
outrA_50=simulate2(L,n,beta,"RA",initial="RRMM",method="RRREG",matlab=matlab,type="Ridge",seed=2017)
outrB_50=simulate2(L,n,beta,"RB",initial="RRMM",method="RRREG",matlab=matlab,type="Ridge",seed=2017)
outrc_50=simulate2(L,n,beta,"RC",initial="RRMM",method="RRREG",matlab=matlab,type="Ridge",seed=2017)
outrD_50=simulate2(L,n,beta,"RD",initial="RRMM",method="RRREG",matlab=matlab,type="Ridge",seed=2017)

#RRREG50(S initial)
outrA_50s=simulate2(L,n,beta,"RA",initial="S",method="RRREG",matlab=matlab,type="Ridge",seed=2017)
outrB_50s=simulate2(L,n,beta,"RB",initial="S",method="RRREG",matlab=matlab,type="Ridge",seed=2017)
outrc_50s=simulate2(L,n,beta,"RC",initial="S",method="RRREG",matlab=matlab,type="Ridge",seed=2017)
outrD_50s=simulate2(L,n,beta,"RD",initial="S",method="RRREG",matlab=matlab,type="Ridge",seed=2017)

#RRREG remove
 

#RRREG0.3
outrA_50.3=simulate2(L,n,beta,"RA",initial="RRMM",method="RRREG",matlab=matlab,type="Ridge",seed=2017,pro=0.3)
outrB_50.3=simulate2(L,n,beta,"RB",initial="RRMM",method="RRREG",matlab=matlab,type="Ridge",seed=2017,pro=0.3)
outrc_50.3=simulate2(L,n,beta,"RC",initial="RRMM",method="RRREG",matlab=matlab,type="Ridge",seed=2017,pro=0.3)
outrD_50.3=simulate2(L,n,beta,"RD",initial="RRMM",method="RRREG",matlab=matlab,type="Ridge",seed=2017,pro=0.3)

matLabDir="D:\\matlab\\RRMM"
matlab=PrepareMatlab(matLabDir)
#------------------------------p=100-------------------------------------#
set.seed(2016)
p=100;pNum=30;p1=20;p2=10
n=50
L=5
u1=runif(pNum,0,1)
u2=runif(pNum,0,1)
sign=ifelse(u2>0.5,1,-1)
beta=c(sign[1:p1]*(2.5+u1[1:p1]),sign[(p1+1):pNum]*u1[(p1+1):pNum],rep(0,p-pNum))


#RRREG_H50
outrA_H50s=simulate2(L,n,beta,"RA",initial="S",method="RRREG",matlab=matlab,type="Ridge",standardize=TRUE,seed=2017)
outrB_H50s=simulate2(L,n,beta,"RB",initial="S",method="RRREG",matlab=matlab,type="Ridge",standardize=TRUE,seed=2017)
outrc_H50s=simulate2(L,n,beta,"RC",initial="S",method="RRREG",matlab=matlab,type="Ridge",standardize=TRUE,seed=2017)
outrD_H50s=simulate2(L,n,beta,"RD",initial="S",method="RRREG",matlab=matlab,type="Ridge",standardize=TRUE,seed=2017)


#RRMM_H50
rA_RRMM_H50=simulate2(L,n,beta,"RA",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017)
rB_RRMM_H50=simulate2(L,n,beta,"RB",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017)
rC_RRMM_H50=simulate2(L,n,beta,"RC",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017)
rD_RRMM_H50=simulate2(L,n,beta,"RD",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2017)
#test
tA_50=simulate2(L,n,beta,"RA",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",seed=2016,
                   lambda1=seq(from=0.2,to=2,length.out=50),lambda2=seq(from=0,to=0,length.out=2))
tB_50=simulate2(L,n,beta,"RB",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",,seed=2016,
                   lambda1=seq(from=0.11,to=5,length.out=100),lambda2=seq(from=0,to=0,length.out=2))
tc_50=simulate2(L,n,beta,"RC",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",,seed=2016,
                   lambda1=seq(from=0.1,to=1,length.out=100),lambda2=seq(from=0,to=0,length.out=2))
tD_50=simulate2(L,n,beta,"RD",method="RRREG",matlab=matlab,useDataFile=FALSE,type="Ridge",,seed=2016,
                   lambda1=seq(from=0.001,to=1,length.out=50),lambda2=seq(from=0,to=0,length.out=2))
set.seed(2016)
out=GenerateDataByModel(n,beta=beta,model="RC",dataType="Ridge")
res=rrreg(x=out$x,y=out$y,initial="RRMM",lambda1=seq(from=0.0001,to=0.001,length.out=30),lambda2=seq(from=0.01,to=2,length.out=30),intercept=TRUE)
GetRobustPe(AddIntercept(out$x),y=out$y,betaHat = res$beta)
res=rrreg(x=out$x,y=out$y,lambda1=seq(from=10,to=20,length.out=100),lambda2=seq(from=0.0001,to=1,length.out=100),intercept="false")

outC=GenerateDataByModel(n,beta=beta,model="A",r=0,dataType="Ridge")

res=rrreg(x=outC$x,y=outC$y,lambda1=1,lambda2=0,intercept="false")

GetRidgeLambda(out$x,out$y,matlab=matlab)

