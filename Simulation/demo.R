
##------------------------------------------------------------------##
##----------------------- Low Dimension ----------------------------##
##------------------------------------------------------------------##
L = 100
n = 50
p = 8
beta = c(3, 2, 1.5, 0, 0, 0, 0, 0)

outA_50 = simulation(L, n, beta, "A", method = "PAWLS", initial = "PAWLS", seed = NULL, useDataFile = TRUE, 
    updateInitial = TRUE)
outB_50 = simulation(L, n, beta, "B", method = "PAWLS", initial = "PAWLS", seed = NULL, useDataFile = TRUE, 
    updateInitial = TRUE)
outC_50 = simulation(L, n, beta, "C", method = "PAWLS", initial = "PAWLS", seed = NULL, useDataFile = TRUE, 
    updateInitial = TRUE)
outD2_50 = simulation(L, n, beta, "D2", method = "PAWLS", initial = "PAWLS", seed = NULL, useDataFile = TRUE, 
    updateInitial = TRUE)

##-------------------------------------------------------------------##
##----------------------- high Dimension ----------------------------##
##-------------------------------------------------------------------##

L = 100
n = 100
p = 500
num = 10
beta = c(rep(2, num), rep(0, p - num))

outA0_500 = simulation(L, n, beta, "A", method = "PAWLS", initial = "PAWLS", seed = 2015, updateInitial = TRUE)
outC0_500 = simulation(L, n, beta, "C", method = "PAWLS", initial = "PAWLS", seed = 2015, updateInitial = TRUE)
outB0_500 = simulation(L, n, beta, "B", method = "PAWLS", initial = "PAWLS", seed = 2015, updateInitial = TRUE)
outD20_500 = simulation(L, n, beta, "D2", method = "PAWLS", initial = "PAWLS", seed = 2015, updateInitial = TRUE)

