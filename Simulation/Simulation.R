simulation = function(L, n, beta = NULL, model = c("A", "B", "C", "D"), p = NULL, method = "PAWLS", 
    matlab = NULL, seed = 2014, useDataFile = FALSE, standardize = FALSE, penalty1 = "1-w0", updateInitial = FALSE, 
    criterion = "BIC", intercept = FALSE, initial = "uniform", range = "cross", type = c("Lasso", 
        "Ridge")) {
    mcount <- length(model)
    
    # define output
    iter = rep(0, L)
    iw = rep(0, L)
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
        ptm = proc.time()
        
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
        b = array(0, dim = c(L, p))
        w = array(0, dim = c(L, n))
        
        # simulation
        for (i in 1:L) {
            # data
            if (useDataFile) {
                out = data[[i]]
            } else {
                out = GenerateDataByModel(n = n, beta = beta, model = model[j], dataType = type)
            }
            
            # try different methods
            if (method == "LAD") {
                init = InitParam(out$x, out$y, method = "LAD")
                b[i, ] = init$beta
            } else if (method == "ROSS") {
                setVariable(matlab, X = out$x)
                setVariable(matlab, y = out$y)
                evaluate(matlab, "[betaRoss]=RossSimulate(X,y)")
                betaRoss = getVariable(matlab, "betaRoss")
                b[i, ] = as.vector(betaRoss$betaRoss)
            } else if (method == "ADL") {
                res = srcdreg(out$x, out$y, penalty1 = penalty1, nlambda1 = 1, initial = "LASSO", 
                  intercept = intercept, standardize = standardize, updateInitialTimes = updateInitialTimes, 
                  criterion = criterion, search = "fixw")
                b[i, ] = res$beta
            } else if (method == "LTS") {
                require(robustHD)
                res = sparseLTS(out$x, out$y, intercept = FALSE)
                b[i, ] = res$coefficients
                w[i, ] = res$wt
            } else if (method == "MMNNG") {
                fbeta = paste("data\\", model[j], n, "X", p, "_beta.rda", sep = "")
                load(fbeta)
                b = betaHat
                break
                # source('mmnngreg.R') res=mmnngreg(out$x,out$y) b[i,]=res$betac[-1]
            } else if (method == "MMNNG_DATA") {
                source("mmnngreg.R")
                res = mmnngreg(out$x, out$y)
                f = paste("data\\", model[j], n, "X", p, ".rda", sep = "")
                fbeta = paste("data\\", model[j], n, "X", p, "_beta.rda", sep = "")
                lf = try(load(f))
                if (class(lf) == "try-error") {
                  data = list(out)
                  betaHat = matrix(nrow = 1, ncol = p, res$betac[-1])
                } else {
                  if (length(data) == L) 
                    stop("data collect finish!")
                  data = c(data, list(out))
                  load(fbeta)
                  betaHat = rbind(betaHat, res$betac[-1])
                }
                save(data, file = f)
                save(betaHat, file = fbeta)
            } else if (method == "PAWLS") {
                if (updateInitial) 
                  updateInitialTimes = 2 else updateInitialTimes = 0
                ptm <- proc.time()
                res = srcdreg(out$x, out$y, penalty1 = penalty1, nlambda1 = 50, nlambda2 = 100, delta = 1e-06, 
                  maxIter = 1000, initial = initial, intercept = intercept, standardize = standardize, 
                  updateInitialTimes = updateInitialTimes, criterion = criterion, search = range)
                times[i] <- (proc.time() - ptm)[1]
                b[i, ] = res$beta
                w[i, ] = res$w
                iter[i] = res$iter
                iw[i] = res$index1
            }
            
            # record result
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
        
        nres[[j]] <- list(model = model[j], CFR = CFR, CFR2 = CFR2, OFR = OFR, PDR = PDR, FDR = FDR, 
            AN = AN, MSE = MSE, TIME = TIME)
    }
    # return
    nres
}
