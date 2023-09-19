
guha_sim1 <- function(nburn=30000,nsamp=20000)
{
  rm(list=ls())
  #if (!require(timeR)) install.packages('timeR')
  
  load("data/test/GuhaData.Rdata") 
  #### Xmat (Predictor Matrix) and y (continuous response vector) are loaded.
  Xmat <- simdata$Xmat
  y <- simdata$y
  ## For later comparisons, load 'true.b'
  true.b <- simdata$true.b
  bin.gen <- simdata$bin.gen
  
  ####  Example values below, and what they mean...
  ####  V <- 20 ## No. of nodes for simulations
  ####  R <- 5 ## Maximum Dimensionality (Varies with simulation settings)
  ####  niter <- 50000  ## No. of iterations
  ####################################################################################
  #### Run the model 
  ####################################################################################
  
  source("additional_materials/test/BNSP-Function.R") #### loads function 'BNSPfunc'
  
  # i added this
  #library(timeR)
  #timer1 <- createTimer()
  #timer1$start("BNSPfunc")
  outlist <- BNSPfunc(Xmat = Xmat, y = y, V = 20, R = 5, niter = nburn+nsamp)
  #timer1$stop("BNSPfunc", comment="outlist <- BNSPfunc (V=20,R=5,niter=50000) finished")
  #getTimer(timer1)
  #### Compute model performance (for inference) metrics.
  ####################################################################################
  ####  ***Post-Analysis Code for Simulated Scenario***      
  ####################################################################################
  
  burnin <- nburn
  postburn <- nsamp
  niter <- nburn+nsamp
  
  #### (1) MSE - Mean Squared Error of edge coefficients
  q <- outlist$q
  betamat <- outlist$betamat
  
  betasample <- matrix(NA,postburn,q)
  for(i in 1:postburn){
    betasample[i,] <- upperTriangle(betamat[i+burnin,,],byrow=T)
  }
  MSE <- mean((colMeans(betasample) - true.b)^2)
  
  
  #### (2) Which Nodes Identified
  rn.gen <- outlist$epsilon_k
  xi_means <- colMeans(rn.gen[(burnin+1):niter,]) #### Detected Nodes if greater than 0.5
  #length(which(colMeans(rn.gen[(burnin+1):niter,])>0.5))
  retlist <- list(gamma=betasample, MSE=MSE, xis=xi_means)
  return(retlist)
}
