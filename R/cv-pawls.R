cv.pawls <-  function(x, y, nfolds=10,seed=123, ...){
  n = length(y)
  p = dim(x)[2]
  if(!is.null(seed)){
    set.seed(seed)
  }
  fit <- pawls(x, y, ...)
  resid_array <- Y <- array(NA, dim=c(n, length(fit$lambda2), length(fit$lambda1)) )
  fold <- sample(1:n %% nfolds)
  fold[fold==0] <- nfolds
  cv.args <- list(...)
  cv.args$lambda1 <- fit$lambda1
  cv.args$lambda2 <- fit$lambda2
  nlambda1 <- length(fit$lambda1)
  nlambda2 <- length(fit$lambda2)
  cl <- makeCluster(parallel::detectCores())
  registerDoParallel(cl)
  res_list <- foreach(i=1:nfolds) %dopar% cvf(i, x, y, fold, cv.args)
  stopCluster(cl)
  registerDoSEQ()
  for (i in 1:nfolds) {
    # res <- cvf(i, x, y, fold, cv.args)
    #resid_array[fold==i,1:res$nlambda2, 1:res$nlambda1] <- res$resid
    # Y[fold==i, 1:res$nlambda2, 1:res$nlambda1] <- res$yhat
    resid_array[fold==i,1:nlambda2, 1:nlambda1] <- res_list[[i]]$resid
    Y[fold==i, 1:nlambda2, 1:nlambda1] <- res_list[[i]]$yhat
  }
  ## Return
  cve <- apply(resid_array, c(2,3), function(r){
    rsq <- r^2
    mean(rsq[rsq<quantile(rsq,probs = 0.8)]) # trimmed mean squared residual
    #mean(rsq)
  })
  opt_loc <- which(cve==min(cve), arr.ind = TRUE)
  
  val <- list(cve=cve, fold=fold, beta=fit$beta[opt_loc[1,1], opt_loc[1,2],],
              w=fit$w[opt_loc[1,1], opt_loc[1,2],],
              lambda1=fit$lambda1, lambda2=fit$lambda2, 
              lambda1.opt=fit$lambda1[opt_loc[1,2]], 
              lambda2.opt=fit$lambda2[opt_loc[1,1]],
              min=min(cve),iter=fit$iter[opt_loc[1,1], opt_loc[1,2]], fit=fit)
  structure(val, class="cv.pawls")
}

cvf <- function(i, XX, y, fold, cv.args) {
  cv.args$x <- XX[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i]
  fit.i <- do.call("pawls", cv.args)
  X2 <- XX[fold==i, , drop=FALSE]
  if(fit.i$intercept){
    X2 <- cbind(rep(1,dim(X2)[1]), X2)
  }
  y2 <- y[fold==i]
  yhat <- apply(fit.i$beta, c(1,2), function(b) X2%*%b)
  resid <- apply(yhat, c(2,3), function(yhat2) y2-yhat2)
  list(resid=resid, nlambda1=length(fit.i$lambda1), nlambda2=length(fit.i$lambda2),
       yhat=yhat)
}

