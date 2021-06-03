bic.pawls <- function(x, y,...){
  res1 <- pawls(x, y, ...)
  res2 <-  BIC_grid(res1$wloss, res1$beta, res1$w)
  fit <- list(beta = res2$beta,
              w = res2$w,
              lambda1 = res1$lambda1,
              lambda2 = res1$lambda2,
              lambda1.opt = res1$lambda1[res2$index1],
              lambda2.opt = res1$lambda2[res2$index2],
              iter = res1$iter,
              ws = res1$w,
              betas = res1$betas,
              raw.bic = res2$raw.bic,
              bic = res2$bic)
  class(fit) <- "bic.pawls"
  fit
}