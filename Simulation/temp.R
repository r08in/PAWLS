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
