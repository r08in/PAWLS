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
