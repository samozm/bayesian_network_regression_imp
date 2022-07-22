
####################################################################################
#### This is the function runs the BNSP model after inputing the data and relevant parameters.
####################################################################################

BNSPfunc <- function(Xmat = Xmat, y = y, V = 20, R = 5, niter = 50000, a.wish = 12){
	
	  #### Load Packages
    pkgs <- c("pscl","LaplacesDemon","MASS","gdata","GIGrvg","mvtnorm","MCMCpack","gtools")
    for (pkg in pkgs)
    {
      if(!require(pkg,character.only = TRUE,quietly = TRUE)) install.packages(pkg)
    }
    
    #### Set Some Values 
    q <- V*(V-1)/2 ## No. of edges
    N <- length(y) ## No. of samples
    a.tau <- 0
    b.tau <- 0
    r <- 1
    delta <- 1 
    a <- 1
    b <- 1
    #a.wish <- 10
    S.scale <- diag(rep(1,R))

    #### Hyperparameters for kappa probability
    a.kap <- c()
    b.kap <- c()
    c.kap <- c()

    for(r in 1:R){
      a.kap[r] <- 1
      b.kap[r] <- r^(1.01)
      c.kap[r] <- 1
    }

     ####################################################################################
     #### Initialize
     ####################################################################################

     betamat <- array(NA,dim=c(niter,V,V))
     mu <- c(rep(NA,niter))
     tau <- c(rep(NA,niter))
     S <- array(NA,dim=c(niter,V,V))
     lambda <- c(rep(NA,niter))
     ulist <- list()
     for (j in 1:niter){
	   ulist[[j]] <- matrix(NA,nrow=V,ncol=R)  ## u1,...,uV for each iteration.
     }
     pi <- numeric()
     xi <- list() ## Same as matrix M
     rn.gen <- matrix(NA,niter,V)
     kappa <- matrix(NA,niter,R)
     pi.kap <- array(NA,dim=c(niter,R,3))

     #### Gibbs Sampler
     #### Initial Values
     tau[1] <- 1
     mu[1] <- 0.8
     for(r in 1:R){
       pi.kap[1,r,] <- c(rdirichlet(1,c(a.kap[r],b.kap[r],c.kap[r])))
       kappa[1,r] <- sample(c(1,0,-1),1,pi.kap[1,r,],replace=T)
     }
     betamat[1,,] <- rnorm(V*V,0,1) 
     S[1,,] <- rgamma(V*V,1,2)
     lambda[1] <- 0.1
     pi[1] <- 0.5
     M <- riwish(a.wish,S.scale)
     xi[[1]] <- M
     #### Each of the 'niter' no. of matrices in the list represents u1,u2,u3,...,uV.(Each of the u's is of dimension R.)
     ulist[[1]] <- matrix(rnorm(V*R,0,1),nrow=V,ncol=R)
     Byatha <- diag(kappa[1,])
     UpU <- ulist[[1]]%*%Byatha%*%t(ulist[[1]])
     #### Upper Triangular portion of U-prime-U
     UpU.up <- upperTriangle(UpU,diag=FALSE,byrow=T)
     mbar <-  UpU.up 
     rn.gen[1,] <- rbinom(V,1,0.5)

     ####################################################################################
     #### Run Model
     ####################################################################################

     for (i in 2:niter){

	    #### 1.Update tau^2
        S.up <- upperTriangle(S[i-1,,],diag=FALSE,byrow=T)
	    betamat.up <- upperTriangle(betamat[i-1,,],diag=FALSE,byrow=T)
	    epsilon <- (y - mu[i-1] -  Xmat%*%betamat.up)
        betamat.u <- betamat.up-mbar
	    tau[i] <- rigamma(1,(a.tau + N/2+V*(V-1)/4),(b.tau + t(epsilon)%*%epsilon/2+sum(betamat.u^2/S.up)/2))
		
	    #### 2.Update uj's (Each u_j for j=1,...,V is a R-dim vector)
        #### Update ulist[[i]].Has V rows for u1,u2,u3,...,uV, and R columns for R-dim uj's.
      
	      mns <- array(NA,dim=c(V,R))
	      cov <- array(NA,dim=c(V,R,R))
        for(j in 1:V){
          if(j==1){
            betabar <- betamat[i-1,1,2:V]
            H <- S[i-1,1,2:V]
          }else if(j==V){
            betabar <- betamat[i-1,1:(V-1),V]
            H <- S[i-1,1:(V-1),V]
          }else{
            betabar <- c(betamat[i-1,1:(j-1),j],betamat[i-1,j,(j+1):V]) 
            H <- c(S[i-1,1:(j-1),j],S[i-1,j,(j+1):V])
          }
          
          ubar <- ulist[[i-1]][-j,]%*%Byatha
          mean.uj <- (solve(t(ubar)%*%diag(1/H)%*%ubar/tau[i] + solve(M)))%*%((t(ubar)%*%diag(1/H)%*%betabar)/tau[i])
    	      cov.uj <- solve(t(ubar)%*%diag(1/H)%*%ubar/tau[i] + solve(M))
    	      pi1 <- pi[i-1]*dmvnorm(betabar,rep(0,length(betabar)),tau[i]*diag(H)+ubar%*%M%*%t(ubar))
    	      pi2 <- (1-pi[i-1])*dmvnorm(betabar,rep(0,length(betabar)),tau[i]*diag(H))
    	      rn.gen[i,j] <- rbinom(1,1,pi1/(pi1+pi2))
    	    if(rn.gen[i,j]==1){
    	      uj <- mvrnorm(1,mean.uj, cov.uj)
    	    }else{uj <- c(rep(0,R))}
          ulist[[i]][j,] <- uj
          mns[j,] <- mean.uj
          cov[j,,] <- cov.uj
        }
    	
      #### 3. Update beta (the parameter 'gamma' in the paper is called 'beta' in this code)
  	  UpU <- ulist[[i]]%*%Byatha%*%t(ulist[[i]])
  	  #### Upper Triangular portion of U-prime-U
  	  UpU.up <- upperTriangle(UpU,diag=FALSE,byrow=T)
  	  mbar <-  UpU.up 	 
  	
  	  Phi.a <- Xmat/sqrt(tau[i])
      alpha.a <- (y-Xmat%*%mbar-mu[i-1])/sqrt(tau[i])
      S.up.a <- tau[i]*S.up
      u.a <- rnorm(length(S.up.a),0,sqrt(S.up.a))
      delta.a <- rnorm(nrow(Xmat))
      v.a <- Phi.a%*%u.a+delta.a
      w.a <- solve(Phi.a%*%diag(S.up.a)%*%t(Phi.a)+diag(nrow(Xmat)))%*%(alpha.a-v.a)
      theta.a <- u.a+diag(S.up.a)%*%t(Phi.a)%*%w.a
      beta.update <- c(theta.a)+mbar
      upperTriangle(betamat[i,,],byrow=T) <- beta.update

      #### 4. Update S
      S.update <- sapply(1:length(beta.update),function(l){rgig(1,1/2,(beta.update[l]-mbar[l])^2/tau[i],lambda[i-1])})
      upperTriangle(S[i,,],byrow=T) <- S.update

      #### 5. Update lambda^2
      lambda[i] <- rgamma(1,r+V*(V-1)/2,sum(S.update)/2+delta)
   
      #### 6. Update pi (corresponding to the U's)
      pi[i] <- rbeta(1,(a+sum(rn.gen[i,])),(b+V-sum(rn.gen[i,])))
    
      #### 7. Update xi
      len.nonzero.ui <- length(which(rn.gen[i,]==1))
      U.U <- matrix(0,R,R)
      for(j in 1:V){ U.U <- U.U+tcrossprod(ulist[[i]][j,])}
      xi[[i]] <- riwish(a.wish+len.nonzero.ui,S.scale+U.U)
      M <- xi[[i]]

      #### 8. Update mu
      mean.mu <- mean(y-Xmat%*%beta.update)
      sd.mu <- sqrt(tau[i]/N)
      mu[i] <- rnorm(1,mean.mu,sd.mu)
    
      #### Byatha means 'Lambda' matrix
      #### 9. Update kappa (The Lambda matrix)
      for(r in 1:R){
        kappa.temp1 <- kappa[i-1,]
        kappa.temp2 <- kappa[i-1,]
        kappa.temp3 <- kappa[i-1,]
        kappa.temp1[r] <- 1
        kappa.temp2[r] <- 0
        kappa.temp3[r] <- -1
        Byth1 <- diag(kappa.temp1)
        Byth2 <- diag(kappa.temp2)
        Byth3 <- diag(kappa.temp3)
        UpU1 <- ulist[[i]]%*%Byth1%*%t(ulist[[i]])
        UpU2 <- ulist[[i]]%*%Byth2%*%t(ulist[[i]])
        UpU3 <- ulist[[i]]%*%Byth3%*%t(ulist[[i]])
        #if(i==32010){
        #  browser()
        #}
  	    ## Upper Triangular portion of U-prime-U
  	    UpU.up1 <- upperTriangle(UpU1,diag=FALSE,byrow=T)
  	    mbar1 <-  UpU.up1 
        UpU.up2 <- upperTriangle(UpU2,diag=FALSE,byrow=T)
	      mbar2 <-  UpU.up2
        UpU.up3 <- upperTriangle(UpU3,diag=FALSE,byrow=T)
	      mbar3 <-  UpU.up3
        prob.up1 <- sum(dnorm(beta.update,mbar1,sqrt(tau[i]*S.update),log=TRUE))
        prob.up2 <- sum(dnorm(beta.update,mbar2,sqrt(tau[i]*S.update),log=TRUE))
        prob.up3 <- sum(dnorm(beta.update,mbar3,sqrt(tau[i]*S.update),log=TRUE))
        prob.f.up <- c(prob.up1,prob.up2,prob.up3)
        ind.prob <- which(prob.f.up==max(prob.f.up))[1]
        if(prob.f.up[ind.prob]==1)
        {
          #print("bot 0")
          #print(prob.f.up)
        }
        puplo <- exp(prob.f.up-prob.f.up[ind.prob])
        
        kappa[i,r] <- sample(c(1,0,-1),1,prob=pi.kap[i-1,r,]*puplo,replace=T)
      
      }
    
      Byatha <- diag(kappa[i,])
      UpU <- ulist[[i]]%*%Byatha%*%t(ulist[[i]])
      #### Upper Triangular portion of U-prime-U
      UpU.up <- upperTriangle(UpU,diag=FALSE,byrow=T)
      mbar <-  UpU.up
    
      #### 10. Update pi.kap
      if(kappa[i,r]==1){
        mm.new <- c(1,0,0)
      }else if(kappa[i,r]==0){
        mm.new <- c(0,1,0)
      }else{mm.new <- c(0,0,1)}
      
      for(r in 1:R){
        pi.kap[i,r,] <- rdirichlet(1,c(a.kap[r],b.kap[r],c.kap[r])+mm.new)
      }
    
      #print(i)
  }
  
  outlist <- list(betamat = betamat, q=q, lambda=kappa, epsilon_k = rn.gen, theta=lambda, Lambda=Byatha,tau=tau, u=ulist, pi=pi.kap, mu=mu, M=xi, delta=pi,S=S)
  return(outlist)
}

