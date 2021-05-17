
guha_sim2_case1 <- function(nburn=30000, nsamp=20000)
{
  rm(list=ls())
  #if (!require(timeR)) install.packages('timeR')
  
  sim2 <- read.csv("data/simulation2_case1.csv")
  sim2_b <- read.csv("data/simulation2_case1_bs.csv")
  
  y <- sim2$y
  Xmat <- data.matrix(sim2[1:190])
  source("test/BNSP-Function.R") #### loads function 'BNSPfunc'
  
  #timer1 <- createTimer()
  #timer1$start("BNSPfunc")
  outlist <- BNSPfunc(Xmat = Xmat, y = y, V = 20, R = 5, niter = nburn+nsamp, a.wish=10)
  #timer1$stop("BNSPfunc", comment="outlist <- BNSPfunc (V=20,R=5,niter=50000) finished")
  #getTimer(timer1)
  #### Compute model performance (for inference) metrics.
  ####################################################################################
  ####  ***Post-Analysis Code for Simulated Scenario***      
  ####################################################################################
  
  true.b <- sim2_b$B
  
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