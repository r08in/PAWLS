source("Simulation/Simulation.R")
# n=100,500
L = 10
n = 100
p = 500
num = 10
beta = c(rep(2, num), rep(0, p - num))


APAWLS_AIC <- simulation(L, n, beta, c("A", "B", "C", "D","E"), method = "PAWLS", initial = "PAWLS",
                          lambda1.min=1e-3, lambda2.min=0.05,
                          seed = 2016, useDataFile = FALSE, updateInitial = FALSE, 
                          intercept = TRUE,initCrit = "AIC" )
APAWLS_BIC <- simulation(L, n, beta, c("A", "B", "C", "D","E"), method = "PAWLS", initial = "PAWLS",
                           lambda1.min=1e-3, lambda2.min=0.05,
                           seed = 2016, useDataFile = FALSE, updateInitial = FALSE, 
                           intercept = TRUE,initCrit = "BIC" )

PAWLS_AIC <- simulation(L, n, beta, c( "A", "B", "C", "D","E"), method = "PAWLS", initial = "uniform",
                   #lambda1.min=1e-3, lambda2.min=0.05,
                   lambda1.min=1e-3, lambda2.min=0.001,
                   seed = 2016, useDataFile = FALSE, updateInitial = FALSE, 
                   intercept = TRUE,criterion = "AIC", search="all")
PAWLS_BIC <- simulation(L, n, beta, c("A", "B", "C", "D","E"), method = "PAWLS", initial = "uniform",
                   #lambda1.min=1e-3, lambda2.min=0.05,
                   lambda1.min=1e-3, lambda2.min=0.05,
                   seed = 2016, useDataFile = FALSE, updateInitial = FALSE, 
                   intercept = TRUE,criterion = "BIC",search="all")

save(APAWLS_AIC , file = "Output/APAWLS_AIC.rda")
save(APAWLS_BIC , file = "Output/APAWLS_BIC.rda")
save(PAWLS_AIC , file = "Output/PAWLS_AIC.rda")
save(PAWLS_BIC , file = "Output/PAWLS_BIC.rda")
load("Output/APAWLs_AIC.rda")
load("Output/APAWLs_BIC.rda")
load("Output/PAWLs_AIC.rda")
load("Output/PAWLs_BIC.rda")

#-----------------------------------------------------------------------

source("Simulation/Simulation.R")
# n=50,p=8----------------------------------------------
L = 100
n = 50
p = 8
beta = c(3, 2, 1.5, 0, 0, 0, 0, 0)
#beta = c(0, 0, 0, 0, 0, 0, 0, 0)


APAWLS_AIC <- simulation(L, n, beta, c("A", "B", "C", "D","E"), method = "PAWLS", initial = "PAWLS",
                         lambda1.min=1e-3, lambda2.min=0.05,
                         seed = 2016, useDataFile = FALSE, updateInitial = FALSE, 
                         intercept = TRUE,initCrit = "AIC" )
APAWLS_BIC <- simulation(L, n, beta, c("A", "B", "C", "D","E"), method = "PAWLS", initial = "PAWLS",
                         lambda1.min=1e-3, lambda2.min=0.05,
                         seed = 2016, useDataFile = FALSE, updateInitial = FALSE, 
                         intercept = TRUE,initCrit = "BIC" )

PAWLS_AIC <- simulation(L, n, beta, c( "A", "B", "C", "D","E"), method = "PAWLS", initial = "uniform",
                        #lambda1.min=1e-3, lambda2.min=0.05,
                        lambda1.min=1e-3, lambda2.min=0.001,
                        seed = 2016, useDataFile = FALSE, updateInitial = FALSE, 
                        intercept = TRUE,criterion = "AIC", search="all")
PAWLS_BIC <- simulation(L, n, beta, c("A", "B", "C", "D","E"), method = "PAWLS", initial = "uniform",
                        #lambda1.min=1e-3, lambda2.min=0.05,
                        lambda1.min=1e-3, lambda2.min=0.05,
                        seed = 2016, useDataFile = FALSE, updateInitial = FALSE, 
                        intercept = TRUE,criterion = "BIC",search="all")

save(APAWLS_AIC , file = "Output/APAWLS_AIC.rda")
save(APAWLS_BIC , file = "Output/APAWLS_BIC.rda")
save(PAWLS_AIC , file = "Output/PAWLS_AIC.rda")
save(PAWLS_BIC , file = "Output/PAWLS_BIC.rda")
load("Output/APAWLs_AIC.rda")
load("Output/APAWLs_BIC.rda")
load("Output/PAWLs_AIC.rda")
load("Output/PAWLs_BIC.rda")

# Find Lambda
lambda1s <- c(0.01, 0.05, 0.001, 0.005,0)
lambda2s <- c(0.01,0.05,0.001,0.005,0)
nlambda1 <- length(lambda1s)
nlambda2 <- length(lambda2s)
overall <- matrix(0,nrow=nlambda1,ncol=nlambda2)
beat <- matrix(0,nrow=nlambda1,ncol=nlambda2)
se <- matrix(0,nrow=nlambda1,ncol=nlambda2)
for(i in 1: length(lambda1s))
  for(j in 1: length(lambda2s)){
    lambda1.min <- lambda1s[i]
    lambda2.min <- lambda2s[j]
    res <- simulation(L, n, beta, c("A", "B", "C", "D","E"), method = "PAWLS", initial = "PAWLS", lambda1.min=lambda1.min, 
                      lambda2.min=lambda2.min, seed = NULL, useDataFile = TRUE, updateInitial = FALSE, intercept = TRUE)
    overall[i,j] <- res[[6]]$overall
    beat[i,j] <- res[[6]]$beat
    se[i,j] <- res[[6]]$se
  }


endreq <- function(res,nlambda,m, t=c(1,2),L=100){
  s=0
  if(t==1){
    for(i in 1:m){
      s <- s+ sum(res[[m]]$iw==nlambda)
    }
  } else{
    for(i in 1:m){
      s <- s+ sum(res[[m]]$ib==nlambda)
    }
  }
  s/(L * m)
}


## generate data
# example is not high-dimensional to keep computation time low
library("mvtnorm")
set.seed(1234)  # for reproducibility
n <- 100  # number of observations
p <- 25   # number of variables
beta <- rep.int(c(1, 0), c(5, p-5))  # coefficients
sigma <- 0.5      # controls signal-to-noise ratio
epsilon <- 0.1    # contamination level
Sigma <- 0.5^t(sapply(1:p, function(i, j) abs(i-j), 1:p))
x <- rmvnorm(n, sigma=Sigma)    # predictor matrix
e <- rnorm(n)                   # error terms
i <- 1:ceiling(epsilon*n)       # observations to be contaminated
e[i] <- e[i] + 5                # vertical outliers
y <- c(x %*% beta + sigma * e)  # response
x[i,] <- x[i,] + 5              # bad leverage points

## fit sparse LTS model for one value of lambda
out_LTS <- sparseLTS(x, y, lambda = 0.05, mode = "fraction")

## fit sparse LTS models over a grid of values for lambda
frac <- seq(0.2, 0.05, by = -0.05)
out_LTS <- sparseLTS(x, y, lambda = frac, mode = "fraction")

