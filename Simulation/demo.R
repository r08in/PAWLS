## generate data
# example is not high-dimensional to keep computation time low
library("mvtnorm")
set.seed(1234)  # for reproducibility
n = 100 # number of observations
p = 8 # number of variables
beta = c(1, 2, 3, 0, 0, 0, 0, 0) # coefficients
sigma <- 1      # controls signal-to-noise ratio
epsilon <- 0.1    # contamination level
Sigma <- 0.5^t(sapply(1:p, function(i, j) abs(i-j), 1:p))
x <- rmvnorm(n, sigma=Sigma)    # predictor matrix
e <- rnorm(n)                   # error terms
i <- 1:ceiling(epsilon*n)       # observations to be contaminated
e[i] <- e[i] +20               # vertical outliers
y <- c(x %*% beta + sigma * e)  # response
#x[i,] <- x[i,] + 5              # bad leverage points

## fit pawls model over a find grid of tuning parameters
out <- pawls(x,y)
adap.cv.pawls <- function(x, y, ...){
  browser()
  init <- cv.pawls(x,y, ... )
  cv.pawls(x, y, beta0=init$beta, w0=init$w,lambda1=init$lambda1,
           lambda2=init$lambda2,...)
}
out <- adap.cv.pawls(x,y)
                     #ambda2=seq(from=5/n, to=1/n, length.out = 10))
out2 <- pawls(x,y,lambda1 = out$opt.lambda1,lambda2 = out$opt.lambda2)
## fit adaptive pawls model over a find grid of tuning parameters
pawls(x,y,lambda1.min = 0.001, lambda2.min = 0.05, initial = "PAWLS")

## fit adaptive pawls model using corss search over a grid of tuning parameters
pawls(x,y,lambda1.min = 0.001, lambda2.min = 0.05, initial = "PAWLS", search="cross")

