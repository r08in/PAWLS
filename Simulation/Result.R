

MySummary = function(beta, weight, pnum, onum) {
    l1 = dim(beta)[1]
    l2 = dim(beta)[2]
    b = matrix(0, nrow = l1, ncol = l2)
    w = matrix(0, nrow = l1, ncol = l2)
    for (i in 1:l1) for (j in 1:l2) {
        b[i, j] = sum(beta[i, j, ] != 0 + 0) + sum(beta[i, j, (1:pnum)] != 0 + 0) * 0.01
        w[i, j] = sum(weight[i, j, ] != 1 + 0) + sum(weight[i, j, (1:onum)] != 1 + 0) * 0.01
    }
    
    list(b = b, w = w)
}

MySummary2 = function(beta, weight, pnum, onum) {
    l1 = length(beta)
    l2 = dim(beta[[1]])[1]
    b = matrix(0, nrow = l1, ncol = l2)
    w = matrix(0, nrow = l1, ncol = l2)
    for (i in 1:l1) for (j in 1:l2) {
        b[i, j] = sum(beta[[i]][j, ] != 0 + 0) + sum(beta[[i]][j, (1:pnum)] != 0 + 0) * 0.01
        w[i, j] = sum(weight[[i]][j, ] != 1 + 0) + sum(weight[[i]][j, 1:onum] != 1 + 0) * 0.01
    }
    
    list(b = b, w = w)
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
        M = sum(temp/num)/n
        JD = sum(temp == 0)/n
    }
    S = sum(apply(w[, (num + 1):m] != 1, 1, sum)/(m - num))/n
    
    list(M = M, S = S, JD = JD)
    
}

## compare SpareLTS and PAWLS
ComparePE = function(x, y, trp = 0.9, ret = 1500) {
    require(robustHD)
    n = length(y)
    # iteration for ret times
    peSumLTS = 0
    peSum = 0
    for (i in 1:ret) {
        # random split data into training set and test set
        trainIndex = sample(1:n, size = n * trp)
        trainx = x[trainIndex, ]
        trainy = y[trainIndex]
        testx = x[-trainIndex, ]
        testy = y[-trainIndex]
        
        # get beta hap and outliers proportion
        resLTS = sparseLTS(trainx, trainy, intercept = TRUE)
        dpLTS = 0.75
        res = srcdreg(trainx, trainy, intercept = TRUE, initial = "LTS")
        dp = sum(res$w == 1)/n
        
        # compute PE for current iteration
        peLTS = (testy - testx %*% resLTS$coefficients[-1] - resLTS$coefficients[1])^2
        pe = (testy - testx %*% res$beta[-1] - res$beta[1])^2
        peSumLTS = peSumLTS + mean(sort(peLTS)[1:(n * (1 - trp) * dpLTS)])
        peSum = peSum + mean(sort(pe)[1:(n * (1 - trp) * dp)])
    }
    
    # mean PE
    list(mpe = peSum/ret, mpeLTS = peSumLTS/ret)
}
