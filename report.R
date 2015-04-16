#####beta MSE#####
#---example1 begin---#
L=100
se1=c(outA_50$se,outA_ADL50$se,outA_LTS50$se,outA_MM50$se,outA_ROSS50$se)
se2=c(outB_50$se,outB_ADL50$se,outB_LTS50$se,outB_MM50$se,outB_ROSS50$se)
se3=c(outC_50$se,outC_ADL50$se,outC_LTS50$se,outC_MM50$se,outC_ROSS50$se)
se4=c(outD_50$se,outD_ADL50$se,outD_LTS50$se,outD_MM50$se,outD_ROSS50$se)

group=c(rep("PAWLS",L),rep("ADL",L),rep("SparseLTS",L),rep("MMNNG",L),rep("SROS",L))
m1=data.frame(se=se1,group=group)
m2=data.frame(se=se2,group=group)
m3=data.frame(se=se3,group=group)
m4=data.frame(se=se4,group=group)
x11()
attach(mtcars)
par(mfrow=c(2,2))
boxplot(se~group,data=m1,main="model A")
boxplot(se~group,data=m2,main="model B")
boxplot(se~group,data=m3,main="model C")
boxplot(se~group,data=m4,main="model D")

#for high dimension
se1=c(outA_500$se,outA_ADL500$se,outA_LTS500$se)
se3=c(outC_500$se,outC_ADL500$se,outC_LTS500$se)
group=c(rep("PWLQ-VS",L),rep("ADL",L),rep("SparseLTS",L))
m1=data.frame(se=se1,group=group)
m3=data.frame(se=se3,group=group)
x11()
attach(mtcars)
par(mfrow=c(1,2))
boxplot(se~group,data=m1,main="model A")
boxplot(se~group,data=m3,main="model C")
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
w1=OutlierSummary(outA_500$w,pro=0)
w2=OutlierSummary(outB_500$w,pro=0)
w3=OutlierSummary(outC_500$w)
W4=OutlierSummary(outD_50$w)
#---example1 end---#

#---example2 begin---#
wLA_PWLQ=OutlierSummary(outLA_PWLQ$w)
wLB_PWLQ=OutlierSummary(outLB_PWLQ$w)
wLC_PWLQ=OutlierSummary(outLC_PWLQ$w)

#---example2 end---#