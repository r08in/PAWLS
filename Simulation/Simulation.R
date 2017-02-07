source("Simulation/CombineData.R")
source("Simulation/SetupMatlab.R")
source("Simulation/GenerateData.R")
source('Simulation/mmnngreg.R')

simulation = function(L, n, beta = NULL, model = c("A", "B", "C", "D"), p = NULL, method = "PAWLS", 
    matlab = NULL, seed = 2014, useDataFile = FALSE, standardize = TRUE, penalty1 = "1-w0", updateInitial = TRUE, 
    criterion = "BIC", intercept = TRUE, initial = "uniform", lambda1.min=1e-03, lambda2.min=0.05, range = "cross", type = c("Lasso", 
        "Ridge")) {
    mcount <- length(model)
    
    # define output
    iter = rep(0, L)
    iw = rep(0, L)
    ib = rep(0,L)
    cfr <- rep(0, L)
    ofr <- rep(0, L)
    cfr2 <- rep(0, L)
    pdr <- rep(0, L)
    fdr <- rep(0, L)
    over_size <- 2
    msize <- rep(0, L)
    mses <- rep(0, L)
    times <- rep(0, L)
    nres <- array(list(), mcount)
    
    pb <- txtProgressBar(1, mcount * L, style = 3)
    for (j in 1:mcount) {
        # for each model
        
        # initiate
        if (!is.null(seed)) {
            set.seed(seed)
        }
        p = length(beta)
        if (useDataFile) {
            f = paste("data\\", model[j], n, "X", p, ".rda", sep = "")
            load(f)
            beta = data[[1]]$beta
            n = length(data[[1]]$y)
        }
        b = array(0, dim = c(L, ifelse(intercept,p+1,p)))
        w = array(0, dim = c(L, n))
        # ROC
        # simulation
        for (i in 1:L) {
            # data
            if (useDataFile) {
                out = data[[i]]
                beta <- out$beta
            } else {
                out = GenerateDataByModel(n = n, beta = beta, model = model[j], dataType = type)
            }
            
            # try different methods
            if (method == "LAD") {
                init = InitParam(out$x, out$y, method = "LAD")
                b[i, ] = init$beta
            } else if (method == "ROSS") {
              if(intercept) out$x <- AddIntercept(out$x)
              setVariable(matlab, X = out$x)
              setVariable(matlab, y = out$y)
              evaluate(matlab, "[betaRoss time]=RossSimulate(X,y)")
              times[i] <- getVariable(matlab, "time")[[1]]
              betaRoss = getVariable(matlab, "betaRoss")
              res <- list(beta=betaRoss$betaRoss)
              b[i, ] = as.vector(betaRoss$betaRoss)
            } else if (method == "ADL") {
              ptm <- proc.time()
              fit1<- ncvreg(out$x, out$y, penalty='lasso')
              beta.init <- coef(fit1)[-1 ,which.min(BIC(fit1))]
              w <- pmin(1e6, 1/abs(beta.init))
              fit2 <- ncvreg(out$x, out$y,penalty='lasso', penalty.factor=w)
              beta.alasso <- coef(fit2)[, which.min(BIC(fit2))]
              times[i] <- (proc.time() - ptm)[1]
              res <- list(beta=beta.alasso)
            } else if (method == "Lasso") {
              ptm <- proc.time()
              fit1<- ncvreg(out$x, out$y, penalty='lasso')
              
              beta.alasso <- coef(fit1)[, which.min(BIC(fit1))]
              times[i] <- (proc.time() - ptm)[1]
              res <- list(beta=beta.alasso)
            }else if (method == "LTS") {
                require(robustHD)
                ptm <- proc.time()
                res = sparseLTS(out$x, out$y, intercept = TRUE)
                times[i] <- (proc.time() - ptm)[1]
                b[i, ] = res$coefficients
                w[i, ] = res$wt
                res$beta <- res$coefficients
            } else if (method == "IPOD") {
              ptm <- proc.time()
              H <- out$x %*% solve(t(out$x)%*%out$x)%*%t(out$x)
              res <- IPOD(out$x,out$y,H, method = "soft")
              times[i] <- (proc.time() - ptm)[1]
              w[i,] <- ifelse(res$gamma == 0, 1, 0)
              res$beta <- rep(0, p + intercept)
            }else if (method == "MMNNG") {
              ptm <- proc.time()
              try_res <- try(mmnngreg(out$x,out$y))
              times[i] <- (proc.time() - ptm)[1]
              if(class(try_res)[1]=="try-error"){
                browser()
                res$beta <- rep(0,ifelse(intercept, p+1, p))
              } else{
                res$beta <- try_res$betac
              }
            } else if (method == "MMNNG_DATA") {
                # load data file and result file
                dfile <- paste("data\\", model[j], n, "X", p, ".rda", sep = "")
                rfile <- paste("data\\", model[j], n, "X", p,"_res", ".rda", sep = "")
                lf <- try(load(dfile))
                if (class(lf) == "try-error") {
                  data <- list(out)
                  res_temp <- list(cfr=rep(0,L), cfr2=rep(0,L),ofr=rep(0,L),pdr=rep(0,L),
                                   fdr=rep(0,L), msize=rep(0,L),mses=rep(0,L),times=rep(0,L),ind=1, count=rep(0,L))
                } else {
                  if (length(data) == L) break
                  data = c(data, list(out))
                  load(rfile)
                  res_temp$ind <- res_temp$ind+1
                }
                save(data, file = dfile)
                save(res_temp, file=rfile)
                
                # fit mmnngreg
                ptm <- proc.time()
                res = mmnngreg(out$x, out$y)
                times[i] <- (proc.time() - ptm)[1]
                res$beta <- res$betac
                
            } else if (method == "PAWLS") {
              updateInitialTimes <- ifelse(updateInitial, 2, 0)
              ptm <- proc.time()
              res = srcdreg(out$x, out$y, penalty1 = penalty1, nlambda1 = 50, nlambda2 = 100, lambda1.min=lambda1.min, lambda2.min=lambda2.min, delta = 1e-06, 
                maxIter = 1000, initial = initial, intercept = intercept, standardize = standardize, 
                updateInitialTimes = updateInitialTimes, criterion = criterion, search = range)
              times[i] <- (proc.time() - ptm)[1]
              b[i, ] = res$beta
              w[i, ] = res$w
              iter[i] = res$iter
              iw[i] = res$index1
              ib[i] = res$index2
            }
            
            # record result
            if(intercept) res$beta <- res$beta[-1]
            true_set <- which(beta != 0)
            pnum <- length(true_set)
            active_set <- which(res$beta != 0)
            msize[i] <- length(active_set)
            common_size <- length(intersect(true_set, active_set))
            cfr[i] <- ifelse(common_size == pnum & msize[i] == pnum, 1, 0)  # correctly fit
            ofr[i] <- ifelse(common_size == pnum & msize[i] > pnum, 1, 0)  # over fit
            cfr2[i] <- ifelse(common_size == pnum & msize[i] <= pnum + over_size, 1, 0)  # over fit by at most 2
            pdr[i] <- common_size/pnum  # positive discover rate
            fdr[i] <- (msize[i] - common_size)/msize[i]  #
            mses[i] <- sum((res$beta - beta)^2)
            setTxtProgressBar(pb, (j - 1) * L + i)
            if(method == "MMNNG_DATA"){
              load(rfile)
              ind <- res_temp$ind
              res_temp$cfr[ind] <- cfr[i]
              res_temp$cfr2[ind] <- cfr2[i]
              res_temp$ofr[ind] <- ofr[i]
              res_temp$pdr[ind] <- pdr[i]
              res_temp$fdr[ind] <- fdr[i]
              res_temp$msize[ind] <- msize[i]
              res_temp$mses[ind] <- mses[i]
              res_temp$times[ind] <- times[i]
              res_temp$count[ind] <- 1
              save(res_temp, file=rfile)
            }
        }
        if(method == "MMNNG_DATA"){
          load(rfile)
          cfr <- res_temp$cfr
          cfr2 <- res_temp$cfr2
          ofr <- res_temp$ofr
          pdr <- res_temp$pdr
          fdr <- res_temp$fdr
          msize <- res_temp$msize
          mses <- res_temp$mses
          times <- res_temp$times
        }
        # compute measurement MSE
        MSE <- round(mean(mses), 3)
        # CFR, OFR, PDR, FDR, AN
        CFR <- round(mean(cfr), 3) * 100
        OFR <- round(mean(ofr), 3) * 100
        CFR2 <- round(mean(cfr2), 3) * 100
        PDR <- round(mean(pdr), 3) * 100
        FDR <- round(mean(fdr), 3) * 100
        AN <- round(mean(msize), 3)
        TIME <- sum(times)
        # outlier dectection
        OD <- "not applicable."
        if(method == "PAWLS" || method == "LTS"|| method=="IPOD" ){
          pro <- 0.1
          if(model[j] == "A" || model[j] == "B")
            pro <- 0
          OD <- OutlierSummary(w, pro)
          roc <- ComputeROC(w,pro=pro)
          OD$tpr <- roc$tpr
          OD$fpr <- roc$fpr
        }
        nres[[j]] <- list(model = model[j], CFR = CFR, CFR2 = CFR2, OFR = OFR, PDR = PDR, FDR = FDR, 
            AN = AN, MSE = MSE, mses=mses, TIME = TIME, iw=iw, ib=ib,OD=OD)
    }
    # return
    # Compute Score
    std.score <- c(74,88,72,63,79)
    overall <- 0
    beat <- 1
    se <- 0
    for(i in 1:length(std.score))
    {
      if(nres[[i]]$CFR-std.score[i]<0) beat <- 0
      overall <- overall+(nres[[i]]$CFR-std.score[i])
      se <- se+(nres[[i]]$CFR-std.score[i])^2
    }
    nres[[mcount+1]] <- list(overall=overall,beat=beat,se=se)
    nres
}

OutlierSummary = function(w, pro = 0.1) {
  n = dim(w)[1]
  m = dim(w)[2]
  num = round(m * pro)
  if (num == 0) {
    M = 0
    JD = 1
  } else {
    temp = apply(w[, 1:num] == 1, 1, sum)
    M = mean(temp/num)
    JD = mean(temp == 0)
  }
  S = mean(apply(w[, (num + 1):m] != 1, 1, sum)/(m - num))
  
  list(M = M, S = S, JD = JD)
  
}

ComputeROC= function(w, cutoff=seq(0,1.01,by=0.01), pro=0.1)
{
  l <- length(cutoff)
  L <- dim(w)[1]
  n <- dim(w)[2]
  ps <- NULL
  onum <- n*pro
  if(onum!=0) ps <-1:onum 
  tpr <- rep(0,l)
  fpr <- rep(0,l)
  for(j in 1 : L){
    for(i in 1 : l){
      ps_temp <- which(w[j,] < cutoff[i])
      onum_temp <- length(ps_temp)
      m <- length(intersect(ps,ps_temp))
      tpr[i] <- tpr[i] + m/ onum 
      fpr[i] <- fpr[i] + (onum_temp - m) / (n - onum)
    }
  }
  tpr <- tpr/L
  fpr <- fpr/L
  list(tpr=tpr,fpr=fpr)
}

