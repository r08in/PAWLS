library(RobustCD)
setwd("D:/RProject/RobustCD")
source("Simulation/GenerateData.R")
source("Simulation/Simulation.R")
source("Simulation/CombineData.R")
source("Simulation/SetupMatlab.R")

context("test simulation result")

test_that("simulation result of APAWLS for Case D ",{
  test.res <-simulation(100,50,c(3,2,1.5,0,0,0,0,0),"D",
                                    method="PAWLS",initial="PAWLS",seed=NULL,
                                    useDataFile=TRUE,updateInitial=TRUE)
  test.res[[1]]$TIME <-0
  expect_equal(test.res[[1]], list(model="D",CFR=70,CFR2=81,OFR=15, 
                                   PDR=95,FDR=8,AN=3.2,MSE=1.383,TIME=0))
})