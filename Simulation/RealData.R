# This file for real data analysis

PrepareData = function(file = "fltr.csv", yname = "1389163_at", screenNum = 500) {
    data = read.table(file, header = TRUE)
    data = t(data)
    y = data[-1, (data[1, ] == yname)]
    x = data[-1, (data[1, ] != yname)]
    # colnames(x)=data[1,(data[1,]!=yname)]
    y = as.numeric(y)
    x = apply(x, 2, as.numeric)
    x = ScreenData(y, x, screenNum)$x
    list(y = y, x = x, data = data)
    
}

ScreenData = function(y, x, screenNum) {
    m = dim(x)[2]
    snum1 = ifelse(m > 4 * screenNum, 4 * screenNum, m)
    index1 = order(apply(x, 2, var), decreasing = TRUE)[1:snum1]
    xx = x[, index1]
    index = order(abs(cor(y, xx)), decreasing = TRUE)[1:screenNum]
    list(y = y, x = xx[, index])
}



paramPlot = function(res, param = "both", label1 = NULL, label2 = NULL, xlimw = NULL, xlimb = NULL, 
    ylimw = NULL, ylimb = NULL) {
    logw = NULL
    loglambdaw = NULL
    b = NULL
    loglambdab = NULL
    nlambdaw = length(res$lambda1s)
    nlambdab = length(res$lambda2s)
    n = length(res$w)
    p = length(res$beta[-1])
    for (i in 1:nlambdaw) {
        logw = c(logw, -log(res$ws[i, ]))
        loglambdaw = c(loglambdaw, rep(-log(res$lambda1s[i]), n))
    }
    
    for (j in 1:(nlambdab - 1)) {
        b = c(b, res$betas[j, -1])
        loglambdab = c(loglambdab, rep(-log(res$lambda2s[j]), p))
    }
    attach(mtcars)
    if (param == "both") {
        # par(mfrow=c(1,2))
        plot(loglambdab, b, xlim = xlimb, ylim = ylimb, xlab = expression(paste("-log(", lambda["1n"], 
            ")")), ylab = expression(hat(beta)[j]), pch = "-")
        # identify(loglambdab,b,labels=label2)
        abline(v = -log(res$lambda2), lty = 2, col = "gray60")
        plot(loglambdaw, logw, xlim = xlimw, ylim = ylimw, xlab = expression(paste("-log(", lambda["2n"], 
            ")")), ylab = expression(-log(hat(w)[i])), pch = "-")
        abline(v = -log(res$lambda1), lty = 2, col = "gray60")
    } else if (param == "1") {
        plot(loglambdaw, logw, xlab = expression(paste("-log(", hat(lambda[2]), ")")), ylab = expression(-log(hat(w))), 
            pch = "-")
        abline(v = -log(res$lambda1), lty = 2)
    } else if (param == "2") {
        plot(loglambdab, b, xlab = expression(paste("-log(", hat(lambda[1]), ")")), ylab = expression(hat(beta)), 
            pch = "-")
        abline(v = -log(res$lambda2), lty = 2)
        # identify(loglambdab,b,labels=label2)
    }
}
