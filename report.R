#####beta MSE#####
#---example1 begin---#
L=100
se1=c(outA_50$se,outA_ADL50$se,outA_LTS50$se,outA_MM50$se,outA_ROSS50$se)
se2=c(outB_50$se,outB_ADL50$se,outB_LTS50$se,outB_MM50$se,outB_ROSS50$se)
se3=c(outC_50$se,outC_LTS50$se,outC_MM50$se,outC_ROSS50$se)
se4=c(outD2_50$se,outD2_LTS50$se,outD2_MM50$se,outD2_ROSS50$se)

group=c(rep("PAWLS",L),rep("ADL",L),rep("SparseLTS",L),rep("MMNNG",L),rep("SROS",L))
group2=c(rep("PAWLS",L),rep("SparseLTS",L),rep("MMNNG",L),rep("SROS",L))
m1=data.frame(se=se1,group=group)
m2=data.frame(se=se2,group=group)
m3=data.frame(se=se3,group=group2)
m4=data.frame(se=se4,group=group2)
x11()
attach(mtcars)
par(mfrow=c(2,2))
boxplot(se~group,data=m1,main="model A",outline=FALSE)
boxplot(se~group,data=m2,main="model B",outline=FALSE)
boxplot(se~group,data=m3,main="model C",outline=FALSE)
boxplot(se~group,data=m4,main="model D",outline=FALSE)

#for high dimension
se1=c(outA_500$se,outA_ADL500$se,outA_LTS500$se)
se2=c(outB_500$se,outB_ADL500$se,outB_LTS500$se)
se3=c(outC_500$se,outC_LTS500$se)
se4=c(outD2_500$se,outD2_LTS500$se)
group=c(rep("PAWL",L),rep("ADL",L),rep("SparseLTS",L))
group2=c(rep("PAWL",L),rep("SparseLTS",L))
m1=data.frame(se=se1,group=group)
m2=data.frame(se=se2,group=group)
m3=data.frame(se=se3,group=group2)
m4=data.frame(se=se4,group=group2)
x11()
attach(mtcars)
par(mfrow=c(2,2))
boxplot(se~group,data=m1,main="model A",outline=FALSE)
boxplot(se~group,data=m2,main="model B",outline=FALSE)
boxplot(se~group,data=m3,main="model C",outline=FALSE)
boxplot(se~group,data=m4,main="model D",outline=FALSE)
#---example1 edn---#

#---example2 begin---#
L=100
se1=c(outLA_PWLQ$se,outLA_MM$se,outLA_LTS$se,outLA_ADL$se)
se2=c(outLB_PWLQ$se,outLB_MM$se,outLB_LTS$se,outLB_ADL$se)
se3=c(outLC_PWLQ$se,outLC_LTS$se,outLC_ADL$se)

group=c(rep("PWLQ",L),rep("MMNNG",L),rep("SpareLTS",L),rep("ADL",L))
m1=data.frame(se=se1,group=group)
m2=data.frame(se=se2,group=group)
group=c(rep("PWLQ",L),rep("Spare LTS",L),rep("ADL",L))
m3=data.frame(se=se3,group=group)
attach(mtcars)
par(mfrow=c(2,2))
boxplot(se~group,data=m1,main="model LA")
boxplot(se~group,data=m2,main="model LB")
boxplot(se~group,data=m3,main="model LC")
#boxplot(se~group,data=m4,main="model D",ylim=c(0,8))
#---example2 end---#
#####beta MSE end#####

##### w #####
#---example1 begin---#
w1=OutlierSummary(outA_50$w,pro=0)
w2=OutlierSummary(outB_50$w,pro=0)
w3=OutlierSummary(outC_50$w)
w4=OutlierSummary(outD2_50$w)

#LTS
wLTS1=OutlierSummary(outA_LTS50$w,pro=0)
wLTS2=OutlierSummary(outB_LTS50$w,pro=0)
wLTS3=OutlierSummary(outC_LTS50$w)
wLTS4=OutlierSummary(outD2_LTS50$w)


wH1=OutlierSummary(outA_500$w,pro=0)
wH2=OutlierSummary(outB_500$w,pro=0)
wH3=OutlierSummary(outC_500$w)
wH4=OutlierSummary(outD2_500$w)

wHLTS1=OutlierSummary(outA_LTS500$w,pro=0)
wHLTS2=OutlierSummary(outB_LTS500$w,pro=0)
wHLTS3=OutlierSummary(outC_LTS500$w)
wHLTS4=OutlierSummary(outD2_LTS50$w)
#---example1 end---#

#---example2 begin---#
wLA_PWLQ=OutlierSummary(outLA_PWLQ$w)
wLB_PWLQ=OutlierSummary(outLB_PWLQ$w)
wLC_PWLQ=OutlierSummary(outLC_PWLQ$w)

#---example2 end---#

##plot air pollution ##

 paramPlot=function(res,param="both")
{
  logw=NULL
  loglambdaw=NULL 
  b=NULL
  loglambdab=NULL
  nlambdaw=length(res$res$lambda1)
  nlambdab=length(res$res$lambda2)
  indexw=res$index1
  indexb=res$index2
  n=length(res$w)
  p=length(res$beta[-1])
  for(i in 1:nlambdaw)
  {
    logw=cbind(logw,-log(res$res$w[i,indexb,]))
    loglambdaw=cbind(loglambdaw,rep(-log(res$res$lambda1[i]),n))
  }
  
  for(j in 1:nlambdab)
  {
    b=cbind(b,res$res$beta[indexw,j,-1])
    loglambdab=cbind(loglambdab,rep(-log(res$res$lambda2[j]),p))
  }
  attach(mtcars)
  if(param=="both")
  {
    par(mfrow=c(1,2))
    plot(loglambdab,b,xlab=expression(paste("-log(",hat(lambda[1]),")")),ylab=expression(hat(beta)),pch="-")
    abline(v=-log(res$lambda2),lty=2)
    plot(loglambdaw,logw,xlab=expression(paste("-log(",hat(lambda[2]),")")),ylab=expression(-log(hat(w))),pch="-")
    abline(v=-log(res$lambda1),lty=2)
  }
  else if(param=="1")
  {
    plot(loglambdaw,logw,xlab=expression(paste("-log(",hat(lambda[2]),")")),ylab=expression(-log(hat(w))),pch="-")
    abline(v=-log(res$lambda1),lty=2)
  }
  else if(param=="2")
  {
    plot(loglambdab,b,xlab=expression(paste("-log(",hat(lambda[1]),")")),ylab=expression(hat(beta)),pch="-")
    abline(v=-log(res$lambda2),lty=2)
  }
}

