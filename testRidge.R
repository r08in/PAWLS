#RRMM
L=100
n=50
p=8
beta=rep(1,p)
#beta=c(3,2,1.5,1.5,1,1,0.5,0.5)
matLabDir="D:\\matlab\\RRMM"
matlab=PrepareMatlab(matLabDir)
outA_RRMM50=simulate2(L,n,beta,"A",method="RRMM",matlab=matlab,useDataFile=FALSE,type="Ridge")


#PWRR
outridgeA_50=simulate2(L,n,beta,"A",method="PWRR",matlab=matlab,useDataFile=FALSE,type="Ridge")
outridgec_50=simulate2(L,n,beta,"C",method="PWRR",matlab=maatlab,useDataFile=FALSE,type="Ridge",
                       lambda1=3.2,lambda2=seq(5,0.01,length.out=100))

out=GenerateDataByModel(n,beta=beta,model="A",dataType="Ridge")
res=rrreg(x=out$x,y=out$y,lambda1=seq(from=0.1,to=1,by=0.1),lambda2=seq(from=0.01,to=1,by=0.1),intercept="false")
res=rrreg(x=out$x,y=out$y,lambda1=seq(from=10,to=20,length.out=100),lambda2=seq(from=0.0001,to=1,length.out=100),intercept="false")

outC=GenerateDataByModel(n,beta=beta,model="A",r=0,dataType="Ridge")

res=rrreg(x=outC$x,y=outC$y,lambda1=1,lambda2=0,intercept="false")
