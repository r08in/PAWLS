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


Sigma=diag(8)
for ( i in 1:8) {
  for (j in 1:i) 
    Sigma[i,j]=Sigma[j,i]=0.5^(i-j)
}

library(MASS)
x=mvrnorm(50,mu=rep(0,8),Sigma=Sigma)


#The leverage point can be observed from the hat matrix
hat= x%*%solve(t(x)%*%x)%*%t(x)
which(diag(hat)>2*sum(diag(hat))/50)  #This will decide where the leverage point is approximately


svdx=svd(x)
u=svdx$u; v=svdx$v; d=svdx$d
## notice that the hat matrix is also equal to u%*%t(u)
round(diag(u%*%t(u)),3)==round(diag(hat),3) 


##Now we can generate newX with some leverage point as follows 
newu=u
newu[1:10,8]=10*u[1:10,8]  # In the direction corresponding to  the smallest eigenvalue of $X$-space by some larger value (say by 1)
#newu[11:15,8]=20*u[11:15,8] 
#newu[21:25,8]=20*u[21:25,8]

#round(t(newu)%*%newu,3) 
newX=newu%*%diag(d)%*%v

#compare newhat and hat; new eigenvalues
newhat=newX%*%solve(t(newX)%*%newX)%*%t(newX)
cbind(svd(newX)$d, d)
cbind( diag(newhat), diag(hat) )


# we can verify that some observations of newX matrix are leverage points by checking

which(diag(newhat)>2*sum(diag(newhat))/50)
which(diag(hat)>2*sum(diag(hat))/50)  #This will decide where the leverage point is approximately