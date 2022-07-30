library(lubridate)
library(dplyr)

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
  max_psrf_xi <- data.frame(psrf=0,pi=0,mu=0,R=0,nu=0,n_microbes=0,samplesize=0,simtype=0)
  max_psrf_gamma <- data.frame(psrf=0,pi=0,mu=0,R=0,nu=0,n_microbes=0,samplesize=0,simtype=0)
  current_max_psrf_xi <- data.frame(psrf=0,pi=0,mu=0,R=0,nu=0,n_microbes=0,samplesize=0,simtype=0)
  current_max_psrf_gamma <- data.frame(psrf=0,pi=0,mu=0,R=0,nu=0,n_microbes=0,samplesize=0,simtype=0)
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
                  if(!file.exists(flnm))
                  {
                    print(flnm)  
                    next
                  }
                  psrf <- read.csv(flnm)
                  if(psrf$max_xi[1] > 1.2)
                  {
                    max_psrf_xi <- rbind(max_psrf_xi,data.frame(psrf=psrf$max_xi[1],pi=pi,mu=mu,R=R[r.i],nu=nu,n_microbes=n_m,samplesize=samplesize,simtype=s))
                  }
                  if(psrf$max_gamma[1] > 1.2)
                  {
                    max_psrf_gamma <- rbind(max_psrf_gamma,data.frame(psrf=psrf$max_gamma[1],pi=pi,mu=mu,R=R[r.i],nu=nu,n_microbes=n_m,samplesize=samplesize,simtype=s))
                  }
                  if(psrf$max_xi[1] > current_max_psrf_xi$psrf[1])
                  {
                    current_max_psrf_xi$psrf <- psrf$max_xi[1]
                    current_max_psrf_xi$pi <- c(pi)
                    current_max_psrf_xi$mu <- c(mu)
                    current_max_psrf_xi$R <- c(R[r.i])
                    current_max_psrf_xi$nu <- c(nu)
                    current_max_psrf_xi$n_microbes <- c(n_m)
                    current_max_psrf_xi$simtype <- c(s)
                    current_max_psrf_xi$samplesize <- c(samplesize)
                  }
                  if(psrf$max_gamma[1] > current_max_psrf_gamma$psrf[1])
                  {
                    current_max_psrf_gamma$psrf <- psrf$max_gamma[1]
                    current_max_psrf_gamma$pi <- c(pi)
                    current_max_psrf_gamma$mu <- c(mu)
                    current_max_psrf_gamma$R <- c(R[r.i])
                    current_max_psrf_gamma$nu <- c(nu)
                    current_max_psrf_gamma$n_microbes <- c(n_m)
                    current_max_psrf_gamma$simtype <- c(s)
                    current_max_psrf_gamma$samplesize <- c(samplesize)
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
                if(psrf$max_xi[1] > 1.2)
                {
                  max_psrf_xi <- rbind(max_psrf_xi,data.frame(psrf=psrf$max_xi[1],pi=pi,mu=mu,R=R[r.i],nu=nu,n_microbes=n_m,samplesize=samplesize,simtype=0))
                }
                if(psrf$max_gamma[1] > 1.2)
                {
                  max_psrf_gamma <- rbind(max_psrf_gamma,data.frame(psrf=psrf$max_gamma[1],pi=pi,mu=mu,R=R[r.i],nu=nu,n_microbes=n_m,samplesize=samplesize,simtype=0))
                }
                if(psrf$max_xi[1] > current_max_psrf_xi$psrf[1])
                {
                  current_max_psrf_xi$psrf <- psrf$max_xi[1]
                  current_max_psrf_xi$pi <- c(pi)
                  current_max_psrf_xi$mu <- c(mu)
                  current_max_psrf_xi$R <- c(R[r.i])
                  current_max_psrf_xi$nu <- c(nu)
                  current_max_psrf_xi$n_microbes <- c(n_m)
                  current_max_psrf_xi$samplesize <- c(samplesize)
                }
                if(psrf$max_gamma[1] > current_max_psrf_gamma$psrf[1])
                {
                  current_max_psrf_gamma$psrf <- psrf$max_gamma[1]
                  current_max_psrf_gamma$pi <- c(pi)
                  current_max_psrf_gamma$mu <- c(mu)
                  current_max_psrf_gamma$R <- c(R[r.i])
                  current_max_psrf_gamma$nu <- c(nu)
                  current_max_psrf_gamma$n_microbes <- c(n_m)
                  current_max_psrf_gamma$samplesize <- c(samplesize)
                }
              }
            }
          }
        }
      }
    }
  }
  return(list(max_psrf_xi=max_psrf_xi,max_psrf_gamma=max_psrf_gamma,worst_xi=current_max_psrf_xi,worst_gamma=current_max_psrf_gamma))
}


#check_psrf(c(1),c("0.0","0.3","0.8"),c("0.8","1.6"),c(5,7,9),c(10),c(8,15,22),
#           "results/simulation/chtc-unrealistic-O0-R7/",c(100,500))

#check_psrf(c(1),c("0.3","0.8"),c("0.8","1.6"),c(7),c(10),c(8,15,22),
#           "results/simulation/chtc-unrealistic-R7-psrf/",c(100,500))

#check_psrf(c(1),c("0.0","0.3","0.8"),c("0.8","1.6"),c(5,7,9),c(10),c(8,15,22),
#           "results/simulation/local/unrealistic-results/",c(100,500))

check_psrf(c(1),c("0.0","0.3","0.8"),c("0.8","1.6"),c(5,7,9),c(10),c(8,15,22),
           "results/simulation/kahan/unrealistic/",c(100,500))

check_psrf(c(1),c("0.0","0.3","0.8"),c("0.8","1.6"),c(5,7,9),c(10),c(8,15,22),
           "results/simulation/wid/unrealistic/",c(100,500))

## 7 only
#check_psrf(c(1),c("0.0","0.3","0.8"),c("0.8","1.6"),c(7),c(10),c(8,15,22),
#           "results/simulation/local/unrealistic-results/",c(100,500))



#check_psrf(2,c(0.3,0.8),c(0.8,1.6),c(9),c(10),c(8,22),
#  "results/simulation/chtc-realistic-O0/",
#  c(500,1000),c("additive_phylo", "additive_random", "interaction_phylo", 
#                "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)

#check_psrf(2,c(0.3,0.8),c(0.8,1.6),c(7),c(10),c(8,22),
#           "results/simulation/chtc-realistic-O0-R7/",
#           c(500,1000),c("additive_phylo", "additive_random", "interaction_phylo", 
#                         "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)

check_psrf(2,c(0.3,0.8),c(0.8,1.6),c(7),c(10),c(8,22),
           "results/simulation/kahan/realistic/",
           c(500,1000),c("additive_phylo", "additive_random", "interaction_phylo", 
                         "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)



all_psrf <- function(simnum,pis,mus,R,nus,n_microbes,inpath,samplesizes,simtypes=c(NULL),pref="",edge_mu=NULL)
{
  tot <- length(pis)*length(mus)*length(R)*length(nus)*length(n_microbes)*length(samplesizes)*length(simtypes)
  psrf_xi <- data.frame(psrf=0,pi=0,mu=0,R=0,nu=0,n_microbes=0,samplesize=0,simtype=0)
  psrf_gamma <- data.frame(psrf=0,pi=0,mu=0,R=0,nu=0,n_microbes=0,samplesize=0,simtype=0)
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
                  if(!file.exists(flnm))
                  {
                    print(flnm)  
                    next
                  }
                  psrf <- read.csv(flnm)
                  
                  psrf_xi <- rbind(psrf_xi,data.frame(psrf=psrf$max_xi[1],
                                           pi=pi,mu=mu,R=R[r.i],
                                           nu=nu,n_microbes=n_m,
                                           samplesize=samplesize,
                                           simtype=s))
                  
                  psrf_gamma <- rbind(psrf_gamma,data.frame(psrf=psrf$max_gamma[1],
                                           pi=pi,mu=mu,R=R[r.i],
                                           nu=nu,n_microbes=n_m,
                                           samplesize=samplesize,
                                           simtype=s))
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
                
                psrf_xi <- rbind(psrf_xi,data.frame(psrf=psrf$max_xi[1],
                                     pi=pi,mu=mu,R=R[r.i],
                                     nu=nu,n_microbes=n_m,
                                     samplesize=samplesize,
                                     simtype=""))
            
                psrf_gamma <- rbind(psrf_gamma,data.frame(psrf=psrf$max_gamma[1],
                                         pi=pi,mu=mu,R=R[r.i],
                                         nu=nu,n_microbes=n_m,
                                         samplesize=samplesize,
                                         simtype=""))
              }
            }
          }
        }
      }
    }
  }
  return(list(gam=psrf_gamma,xi=psrf_xi))
}


a_p <- all_psrf(2,c(0.3,0.8),c(0.8,1.6),c(7),c(10),c(8,22),
           "results/simulation/chtc-realistic-O0-R7/",
           c(500,1000),c("additive_phylo", "additive_random", "interaction_phylo", 
                         "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)


psrf_xi <- a_p$xi
psrf_xi %>% filter(psrf_xi$simtype=="additive_phylo" & psrf_xi$samplesize==500) %>% filter(psrf == max(psrf))
psrf_xi %>% filter(psrf_xi$simtype=="additive_random" & psrf_xi$samplesize==500) %>% filter(psrf == max(psrf))
psrf_xi %>% filter(psrf_xi$simtype=="interaction_phylo" & psrf_xi$samplesize==500) %>% filter(psrf == max(psrf))
psrf_xi %>% filter(psrf_xi$simtype=="interaction_random" & psrf_xi$samplesize==500) %>% filter(psrf == max(psrf))
psrf_xi %>% filter(psrf_xi$simtype=="redundant_phylo" & psrf_xi$samplesize==500) %>% filter(psrf == max(psrf))
psrf_xi %>% filter(psrf_xi$simtype=="redundant_random" & psrf_xi$samplesize==500) %>% filter(psrf == max(psrf))

psrf_xi %>% filter(psrf_xi$simtype=="additive_phylo" & psrf_xi$samplesize==1000) %>% filter(psrf == max(psrf))
psrf_xi %>% filter(psrf_xi$simtype=="additive_random" & psrf_xi$samplesize==1000) %>% filter(psrf == max(psrf))
psrf_xi %>% filter(psrf_xi$simtype=="interaction_phylo" & psrf_xi$samplesize==1000) %>% filter(psrf == max(psrf))
psrf_xi %>% filter(psrf_xi$simtype=="interaction_random" & psrf_xi$samplesize==1000) %>% filter(psrf == max(psrf))
psrf_xi %>% filter(psrf_xi$simtype=="redundant_phylo" & psrf_xi$samplesize==1000) %>% filter(psrf == max(psrf))
psrf_xi %>% filter(psrf_xi$simtype=="redundant_random" & psrf_xi$samplesize==1000) %>% filter(psrf == max(psrf))


psrf_gamma <- a_p$gam
psrf_gamma %>% filter(simtype=="additive_phylo" & samplesize==500) %>% filter(psrf == max(psrf))
psrf_gamma %>% filter(simtype=="additive_random" & samplesize==500) %>% filter(psrf == max(psrf))
psrf_gamma %>% filter(simtype=="interaction_phylo" & samplesize==500) %>% filter(psrf == max(psrf))
psrf_gamma %>% filter(simtype=="interaction_random" & samplesize==500) %>% filter(psrf == max(psrf))
psrf_gamma %>% filter(simtype=="redundant_phylo" & samplesize==500) %>% filter(psrf == max(psrf))
psrf_gamma %>% filter(simtype=="redundant_random" & samplesize==500) %>% filter(psrf == max(psrf))

psrf_gamma %>% filter(simtype=="additive_phylo" & samplesize==1000) %>% filter(psrf == max(psrf))
psrf_gamma %>% filter(simtype=="additive_random" & samplesize==1000) %>% filter(psrf == max(psrf))
psrf_gamma %>% filter(simtype=="interaction_phylo" & samplesize==1000) %>% filter(psrf == max(psrf))
psrf_gamma %>% filter(simtype=="interaction_random" & samplesize==1000) %>% filter(psrf == max(psrf))
psrf_gamma %>% filter(simtype=="redundant_phylo" & samplesize==1000) %>% filter(psrf == max(psrf))
psrf_gamma %>% filter(simtype=="redundant_random" & samplesize==1000) %>% filter(psrf == max(psrf))


a_p_unreal <- all_psrf(1,c(0.3,0.8),c(0.8,1.6),c(7),c(10),c(8,15,22),
                "results/simulation/wid/unrealistic/",
                c(500,100))


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


#check_fl_exist(c(1),c("0.0","0.3","0.8"),c("0.8","1.6"),c(5,9),c(10),c(8,15,22),"results/simulation/chtc-unrealistic/",c(100,500))


check_fl_exist(c(1),c("0.0","0.3","0.8"),c("0.8","1.6"),c(5,7,9),c(10),c(8,15,22),
           "results/simulation/chtc-unrealistic-O0-R7/",c(100,500))

#check_fl_exist(2,c(0.3,0.8),c(0.8,1.6),c(9),c(10),c(8,22),
#             "results/simulation/chtc-realistic/",
#             c(500,1000),c("additive_phylo", "additive_random", "interaction_phylo", 
#                          "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)

check_fl_exist(2,c(0.3,0.8),c(0.8,1.6),c(7),c(10),c(8,22),
               "results/simulation/chtc-realistic-O0-R7/",
               c(500,1000),c("additive_phylo", "additive_random", "interaction_phylo", 
                             "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)

check_fl_exist(2,c(0.3,0.8),c(0.8,1.6),c(7),c(10),c(8,22),
               "results/simulation/kahan/chtc-realistic-O0-R7/",
               c(500,1000),c("additive_phylo", "additive_random", "interaction_phylo", 
                             "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)



check_time <- function(simnum,pis,mus,R,nus,n_microbes,inpath,samplesizes,simtypes=c(NULL),pref="",edge_mu=NULL)
{
  tot <- length(pis)*length(mus)*length(R)*length(nus)*length(n_microbes)*length(samplesizes)*length(simtypes)
  tmz <- data.frame(time=0,pi=0,mu=0,R=0,nu=0,n_microbes=0,samplesize=0,simtype=0)
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
                  flnm <- make.filename(inpath,simnum,pi,mu,R[r.i],nu,n_m,"time",samplesize,s,edge_mu) 
                  if(!file.exists(flnm)){next}
                  tm <- read.csv(flnm)
                  tmz <- rbind(tmz,data.frame(time=tm$time,
                                              pi=pi,mu=mu,R=R[r.i],
                                              nu=nu,n_microbes=n_m,
                                              samplesize=samplesize,
                                              simtype=s))
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
                flnm <- make.filename(inpath,simnum,pi,mu,R[r.i],nu,n_m,"time",samplesize) 
                if(!file.exists(flnm)){next}
                tm <- read.csv(flnm)
                tmz <- rbind(tmz,data.frame(time=tm$time,pi=pi,mu=mu,R=R[r.i],
                                            nu=nu,n_microbes=n_m,
                                            samplesize=samplesize,simtype=0))
              }
            }
          }
        }
      }
    }
  }
  return(tmz)
}
tms_u <- check_time(c(1),c("0.0","0.3","0.8"),c("0.8","1.6"),c(5,7,9),c(10),c(8,15,22),
           "results/simulation/local/unrealistic-results/",c(100,500))

seconds_to_period(mean(tms_u[tms_u$samplesize == 100 & 
                            !(tms_u$mu == 1.6 & tms_u$pi == 0.0 & tms_u$n_microbes == 22 & tms_u$R == 9) &
                            !(tms_u$mu == 1.6 & tms_u$pi == 0.3 & tms_u$n_microbes == 22 & tms_u$R == 5) &
                            !(tms_u$mu == 1.6 & tms_u$pi == 0.3 & tms_u$n_microbes == 22 & tms_u$R == 7) &
                            !(tms_u$mu == 1.6 & tms_u$pi == 0.3 & tms_u$n_microbes == 22 & tms_u$R == 9) &
                            !(tms_u$mu == 1.6 & tms_u$pi == 0.8 & tms_u$n_microbes == 22 & tms_u$R == 5) ,
                            "time"]))

seconds_to_period(mean(tms_u[tms_u$samplesize == 500 & 
                            !(tms_u$mu == 0.8 & tms_u$pi == 0.3 & tms_u$n_microbes == 15 & tms_u$R==5) &
                            !(tms_u$mu == 0.8 & tms_u$pi == 0.3 & tms_u$n_microbes == 15 & tms_u$R==7) &
                            !(tms_u$mu == 1.6 & tms_u$pi == 0.3 & tms_u$n_microbes == 8  & tms_u$R==5) ,
                             "time"]))

#seconds_to_period(mean(tms_u[tms_u$samplesize == 100 & tms_u$R == 5 &
#                               !((tms_u$mu==1.6 & tms_u$pi==0.3 & tms_u$n_microbes==22 &
#                                    tms_u$samplesize==100 & tms_u$R==5)), "time"])/2)
#seconds_to_period(mean(tms_u[tms_u$samplesize == 100 & tms_u$R == 7 & 
#                               !((tms_u$mu==1.6 & tms_u$pi==0.3 & tms_u$n_microbes==22 &
#                                    tms_u$samplesize==100 & tms_u$R==7)), "time"])/2)
#seconds_to_period(mean(tms_u[tms_u$samplesize == 100 & tms_u$R == 9 & 
#                               !((tms_u$mu==1.6 & tms_u$pi==0.8 & tms_u$n_microbes==22 &
#                                    tms_u$samplesize==100 & tms_u$R==9)), "time"])/2)
#seconds_to_period(mean(tms_u[tms_u$samplesize == 500 & tms_u$R == 5 & 
#                             !((tms_u$mu==0.8 & tms_u$pi==0.3 & tms_u$n_microbes==15 &
#                                  tms_u$samplesize==500 & tms_u$R==5) | 
#                              (tms_u$mu==0.8 & tms_u$pi==0.3 & tms_u$n_microbes==22 &
#                                  tms_u$samplesize==100 & tms_u$R==5)),"time"])/2)
#seconds_to_period(mean(tms_u[tms_u$samplesize == 500 & tms_u$R == 7, "time"])/2)
#seconds_to_period(mean(tms_u[tms_u$samplesize == 500 & tms_u$R == 9, "time"])/2)



tms <- check_time(2,c(0.3,0.8),c(0.8,1.6),c(7),c(10),c(8,22),
           "results/simulation/chtc-realistic-O0-R7/",
           c(500,1000),c("additive_phylo", "additive_random", "interaction_phylo", 
                         "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)

seconds_to_period(mean(tms[tms$samplesize == 500 & 
                           !((tms$n_microbes==22 & tms$pi==0.8 & 
                               tms$mu==1.6 & tms$samplesize==500 &
                               tms$simtype == "redundant_random") | 
                               (tms$n_microbes==22 & tms$pi==0.8 & 
                                  tms$mu==1.6 & tms$samplesize==500 &
                                  tms$simtype == "redundant_phylo") | 
                               (tms$n_microbes==22 & tms$pi==0.3 & 
                                tms$mu==0.8 & tms$samplesize==500 &
                                tms$simtype == "redundant_phylo"))  & 
                            tms$samplesize != 0, "time"]))
seconds_to_period(mean(tms[tms$samplesize == 1000 & 
                           !((tms$n_microbes==8 & tms$pi==0.3 & 
                               tms$mu==0.8 & tms$samplesize==1000 &
                               tms$simtype == "additive_random") | 
                               (tms$n_microbes==8 & tms$pi==0.3 & 
                                  tms$mu==0.8 & tms$samplesize==1000 &
                                  tms$simtype == "additive_phylo") |
                               (tms$n_microbes==22 & tms$pi==0.3 & 
                                tms$mu==0.8 & tms$samplesize==1000 &
                                tms$simtype == "redundant_phylo")) & 
                           tms$samplesize != 0, "time"]))




mu<-1.6
pi<-0.3
n_microbes<-22
flnm<-sprintf("data/simulation/realistic/edge_mu=0.4_mu=%s_n_microbes=%s_out=XYs_pi=%s_samplesize=500_simnum=2_type=interaction_random.csv",
        mu,n_microbes,pi)
a1.9 <- read.csv(flnm)

n_microbes<-8
flnm<-sprintf("data/simulation/realistic/edge_mu=0.4_mu=%s_n_microbes=%s_out=XYs_pi=%s_samplesize=500_simnum=2_type=interaction_random.csv",
              mu,n_microbes,pi)
b1.9 <- read.csv(flnm)


ggplot() + geom_histogram(alpha=0.5,data=a1.9,aes(x=y,fill="red")) + 
  geom_histogram(alpha=0.5,data=b1.9,aes(x=y,fill="blue")) +
  geom_vline(xintercept=7)


mu<-1.6
pi<-0.8
n_microbes<-22
flnm<-sprintf("data/simulation/realistic/edge_mu=0.4_mu=%s_n_microbes=%s_out=XYs_pi=%s_samplesize=500_simnum=2_type=interaction_random.csv",
              mu,n_microbes,pi)
a2.4 <- read.csv(flnm)

n_microbes<-8
flnm<-sprintf("data/simulation/realistic/edge_mu=0.4_mu=%s_n_microbes=%s_out=XYs_pi=%s_samplesize=500_simnum=2_type=interaction_random.csv",
              mu,n_microbes,pi)
b2.4 <- read.csv(flnm)


ggplot() + geom_histogram(alpha=0.5,data=a2.4,aes(x=y,fill="red")) + 
  geom_histogram(alpha=0.5,data=b2.4,aes(x=y,fill="blue")) +
  geom_vline(xintercept=30)

mu<-0.8
pi<-0.3
n_microbes<-22
flnm<-sprintf("data/simulation/realistic/edge_mu=0.4_mu=%s_n_microbes=%s_out=XYs_pi=%s_samplesize=500_simnum=2_type=interaction_random.csv",
              mu,n_microbes,pi)
a1.1 <- read.csv(flnm)

n_microbes<-8
flnm<-sprintf("data/simulation/realistic/edge_mu=0.4_mu=%s_n_microbes=%s_out=XYs_pi=%s_samplesize=500_simnum=2_type=interaction_random.csv",
              mu,n_microbes,pi)
b1.1 <- read.csv(flnm)


ggplot() + geom_histogram(alpha=0.5,data=a1.1,aes(x=y,fill="red")) + 
  geom_histogram(alpha=0.5,data=b1.1,aes(x=y,fill="blue")) +
  geom_vline(xintercept=3)


mu<-0.8
pi<-0.8
n_microbes<-22
flnm<-sprintf("data/simulation/realistic/edge_mu=0.4_mu=%s_n_microbes=%s_out=XYs_pi=%s_samplesize=500_simnum=2_type=interaction_random.csv",
              mu,n_microbes,pi)
a1.6 <- read.csv(flnm)

n_microbes<-8
flnm<-sprintf("data/simulation/realistic/edge_mu=0.4_mu=%s_n_microbes=%s_out=XYs_pi=%s_samplesize=500_simnum=2_type=interaction_random.csv",
              mu,n_microbes,pi)
b1.6 <- read.csv(flnm)


ggplot() + geom_histogram(alpha=0.5,data=a1.6,aes(x=y,fill="red")) + 
  geom_histogram(alpha=0.5,data=b1.6,aes(x=y,fill="blue")) +
  geom_vline(xintercept=22)



