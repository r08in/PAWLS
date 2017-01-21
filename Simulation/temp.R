L = 100
n = 50
p = 8
beta = c(3, 2, 1.5, 0, 0, 0, 0, 0)
out <- GenerateDataByModel(n = n, beta = beta, model = "C")

install.packages("leapp")
H <- out$x %*% solve(t(out$x)%*%out$x)%*%t(out$x)
IPOD(out$x,out$y,H, method = "soft")

source("https://bioconductor.org/biocLite.R")
biocLite("sva")

par(mfrow=c(1,2))
plot(density(xx[,1]), ylim=c(0,1)); lines(density(out$x[,1]), col="red")
plot(density(xx[,2]), ylim=c(0,1)); lines(density(out$x[,2]), col="red")
x11();
par(mfrow=c(1,2))
plot(density(xx[,3]), ylim=c(0,1)); lines(density(out$x[,3]), col="red")
plot(density(xx[,4]), ylim=c(0,1)); lines(density(out$x[,4]), col="red")
x11();
par(mfrow=c(1,2))
plot(density(xx[,5]), ylim=c(0,1)); lines(density(out$x[,5]), col="red")
plot(density(xx[,6]), ylim=c(0,1)); lines(density(out$x[,6]), col="red")
x11();
par(mfrow=c(1,2))
plot(density(xx[,7]), ylim=c(0,1)); lines(density(out$x[,7]), col="red")
plot(density(xx[,8]), ylim=c(0,1)); lines(density(out$x[,8]), col="red")