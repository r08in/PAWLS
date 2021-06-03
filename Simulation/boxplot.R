Lres <- list(Lres_ADL, Lres_mmnng, Lres_PAMLS, Lres_APAMLS, Lres_LTS, Lres_PAWLS, Lres_APAWLS)
mse_array <- array(NA, dim=c(7, 5, 100))
for(i in 1:length(Lres)){
  res <- Lres[[i]]
  temp <- t(sapply(res, getElement, "mses"))
  mse_array[i,,] <- temp
}

L_method <- c("MMNNG", "SROS-2", "ASROS-2", "SLTS", "PAWLS", "APAWLS")
cases <- c("A", "B", "C", "D", "E")
par(mfrow=c(1,1))
boxplot(t(mse_array[,1,]), outline=FALSE, names=c("ALASSO",L_method), ylab="MSE", main=cases[1])
par(mfrow=c(2,2))
for(i in 2:5){
  boxplot(t(mse_array[2:7,i,]), outline=FALSE, names=L_method, ylab="MSE", main=cases[i])
}


Hres <- list(Hres_ADL, Hres_PAMLS, Hres_APAMLS, Hres_LTS, Hres_PAWLS, Hres_APAWLS)
mse_array <- array(NA, dim=c(6, 5, 100))
for(i in 1:length(Hres)){
  res <- Hres[[i]]
  temp <- t(sapply(res, getElement, "mses"))
  mse_array[i,,] <- temp
}

H_method <- c( "SROS-2", "ASROS-2", "SLTS", "PAWLS", "APAWLS")
cases <- c("A", "B", "C", "D", "E")
par(mfrow=c(1,1))
boxplot(t(mse_array[,1,]), outline=FALSE, names=c("ALASSO",H_method), ylab="MSE", main=cases[1])
par(mfrow=c(2,2))
for(i in 2:5){
  boxplot(t(mse_array[2:6,i,]), outline=FALSE, names=H_method, ylab="MSE", main=cases[i])
}