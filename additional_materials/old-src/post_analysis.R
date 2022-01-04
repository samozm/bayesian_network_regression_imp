make.filename <- function(path,simnum,pi,mu,R,nu,n_microbes,type,samplesize,simtype=NULL,edge_mu=NULL)
{
  if(is.null(simtype))
  {
    return(sprintf("%sR=%s_mu=%s_n_microbes=%s_nu=%s_out=%s_pi=%s_samplesize=%s_simnum=%s.csv",
                   path,R,mu,n_microbes,nu,type,pi,samplesize,simnum))
  }
  return(sprintf("%sR=%s_edge_mu=%s_mu=%s_n_microbes=%s_nu=%s_out=%s_pi=%s_samplesize=%s_simnum=%s_type=%s.csv",
                 path,R,edge_mu,mu,n_microbes,nu,type,pi,samplesize,simnum,simtype))
  
}

check_psrf <- function(simnum,pis,mus,R,nus,n_microbes,inpath,samplesizes,simtypes=c(NULL),pref="",edge_mu=NULL)
{
  tot <- length(pis)*length(mus)*length(R)*length(nus)*length(n_microbes)*length(samplesizes)*length(simtypes)
  mean_psrf_xi <- data.frame(psrf=0,pi=0,mu=0,R=0,nu=0,n_microbes=0,samplesize=0,simtype=0)
  max_psrf_xi <- data.frame(psrf=0,pi=0,mu=0,R=0,nu=0,n_microbes=0,samplesize=0,simtype=0)
  mean_psrf_gamma <- data.frame(psrf=0,pi=0,mu=0,R=0,nu=0,n_microbes=0,samplesize=0,simtype=0)
  max_psrf_gamma <- data.frame(psrf=0,pi=0,mu=0,R=0,nu=0,n_microbes=0,samplesize=0,simtype=0)
  if(!is.null(simtypes)){
    for(s in simtypes)
    { 
      for(samplesize in samplesizes)
      {
        for(r.i in 1:length(R))
        {
          for(n_m in n_microbes)
          {
            for(pi in pis)
            {
              for(mu in mus)
              {
                for(nu in nus)
                {
                  flnm <- make.filename(inpath,simnum,pi,mu,R[r.i],nu,n_m,"psrf",samplesize,s,edge_mu) 
                  if(!file.exists(flnm)){next}
                  psrf <- read.csv(flnm)
                  if(psrf$mean_xi[1] > 1.2)
                  {
                    mean_psrf_xi <- rbind(mean_psrf_xi,data.frame(psrf=psrf$mean_xi[1],pi=pi,mu=mu,R=R[r.i],nu=nu,n_microbes=n_m,samplesize=samplesize,simtype=s))
                  }
                  if(psrf$mean_gamma[1] > 1.2)
                  {
                    mean_psrf_gamma <- rbind(mean_psrf_gamma,data.frame(psrf=psrf$mean_gamma[1],pi=pi,mu=mu,R=R[r.i],nu=nu,n_microbes=n_m,samplesize=samplesize,simtype=s))
                  }
                  if(psrf$max_xi[1] > 1.2)
                  {
                    max_psrf_xi <- rbind(max_psrf_xi,data.frame(psrf=psrf$max_xi[1],pi=pi,mu=mu,R=R[r.i],nu=nu,n_microbes=n_m,samplesize=samplesize,simtype=s))
                  }
                  if(psrf$max_gamma[1] > 1.2)
                  {
                    max_psrf_gamma <- rbind(max_psrf_gamma,data.frame(psrf=psrf$max_gamma[1],pi=pi,mu=mu,R=R[r.i],nu=nu,n_microbes=n_m,samplesize=samplesize,simtype=s))
                  }
                }
              }
            }
          }
        }
      }
    }
  }else{
    for(samplesize in samplesizes)
    {
      for(r.i in 1:length(R))
      {
        for(n_m in n_microbes)
        {
          for(pi in pis)
          {
            for(mu in mus)
            {
              for(nu in nus)
              {
                flnm <- make.filename(inpath,simnum,pi,mu,R[r.i],nu,n_m,"psrf",samplesize) 
                if(!file.exists(flnm)){next}
                psrf <- read.csv(flnm)
                if(psrf$mean_xi[1] > 1.2)
                {
                  mean_psrf_xi <- rbind(mean_psrf_xi,data.frame(psrf=psrf$mean_xi[1],pi=pi,mu=mu,R=R[r.i],nu=nu,n_microbes=n_m,samplesize=samplesize,simtype=0))
                }
                if(psrf$mean_gamma[1] > 1.2)
                {
                  mean_psrf_gamma <- rbind(mean_psrf_gamma,data.frame(psrf=psrf$mean_gamma[1],pi=pi,mu=mu,R=R[r.i],nu=nu,n_microbes=n_m,samplesize=samplesize,simtype=0))
                }
                if(psrf$max_xi[1] > 1.2)
                {
                  max_psrf_xi <- rbind(max_psrf_xi,data.frame(psrf=psrf$max_xi[1],pi=pi,mu=mu,R=R[r.i],nu=nu,n_microbes=n_m,samplesize=samplesize,simtype=0))
                }
                if(psrf$max_gamma[1] > 1.2)
                {
                  max_psrf_gamma <- rbind(max_psrf_gamma,data.frame(psrf=psrf$max_gamma[1],pi=pi,mu=mu,R=R[r.i],nu=nu,n_microbes=n_m,samplesize=samplesize,simtype=0))
                }
              }
            }
          }
        }
      }
    }
  }
  return(list(mean_psrf_xi=mean_psrf_xi,mean_psrf_gamma=mean_psrf_gamma,max_psrf_xi=max_psrf_xi,max_psrf_gamma=max_psrf_gamma))
}


check_psrf(c(1),c("0.3","0.8"),c("0.8","1.6"),c(5,9),c(10),c(8,15,22),"results/simulation/chtc-unrealistic/",c(100,500))



check_psrf(2,c(0.3,0.8),c(0.8,1.6),c(9),c(10),c(8,22),
  "results/simulation/chtc-realistic/",
  c(500,1000),c("additive_phylo", "additive_random", "interaction_phylo", 
                "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)





check_fl_exist <- function(simnum,pis,mus,R,nus,n_microbes,inpath,samplesizes,simtypes=c(NULL),pref="",edge_mu=NULL)
{
  ct <- 0
  if(!is.null(simtypes)){
    for(s in simtypes)
    { 
      for(samplesize in samplesizes)
      {
        for(r.i in 1:length(R))
        {
          for(n_m in n_microbes)
          {
            for(pi in pis)
            {
              for(mu in mus)
              {
                for(nu in nus)
                {
                  flnm <- make.filename(inpath,simnum,pi,mu,R[r.i],nu,n_m,"psrf",samplesize,s,edge_mu) 
                  if(!file.exists(flnm)){
                    print(flnm)
                    ct <- ct + 1
                  }
                }
              }
            }
          }
        }
      }
    }
  }else{
    for(samplesize in samplesizes)
    {
      for(r.i in 1:length(R))
      {
        for(n_m in n_microbes)
        {
          for(pi in pis)
          {
            if(R[r.i] == 5 & pi == "0.0"){next}
            for(mu in mus)
            {
              for(nu in nus)
              {
                flnm <- make.filename(inpath,simnum,pi,mu,R[r.i],nu,n_m,"psrf",samplesize) 
                if(!file.exists(flnm)){
                  print(flnm)
                  ct <- ct+1
                }
              }
            }
          }
        }
      }
    }
  }
  print(ct)
}


check_fl_exist(c(1),c("0.0","0.3","0.8"),c("0.8","1.6"),c(5,9),c(10),c(8,15,22),"results/simulation/chtc-unrealistic/",c(100,500))


check_fl_exist(2,c(0.3,0.8),c(0.8,1.6),c(9),c(10),c(8,22),
             "results/simulation/chtc-realistic/",
             c(500,1000),c("additive_phylo", "additive_random", "interaction_phylo", 
                          "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)





