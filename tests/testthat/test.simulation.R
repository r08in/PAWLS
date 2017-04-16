library(RobustCD)
#setwd("E:/R/RobustCD")
setwd("D:/RProject/RobustCD")
source("Simulation/GenerateData.R")
source("Simulation/Simulation.R")
source("Simulation/CombineData.R")
source("Simulation/SetupMatlab.R")

context("test simulation result")

test_that("simulation result of APAWLS for Case D ",{
  test.res <-simulation(20, 50, c(3, 2, 1.5, 0, 0, 0, 0, 0), c("D"), method = "PAWLS", initial = "PAWLS", 
                        lambda1.min=0.05, lambda2.min=0.001,
                        seed = NULL, useDataFile = TRUE,       
                        updateInitial =FALSE, intercept = TRUE)
  expect_equal(test.res[[1]][1:5], list(model="D",CFR=90,CFR2=90,OFR=0, 
                                   PDR=96.7))
})