L = 100
n = 50
p = 8
beta = c(3, 2, 1.5, 0, 0, 0, 0, 0)
out <- GenerateDataByModel(n = n, beta = beta, model = "C")

install.packages("leapp")
H <- out$x %*% solve(t(out$x)%*%out$x)%*%t(out$x)
IPOD(out$x,out$y,H, method = "soft")
