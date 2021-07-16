
library(ggplot2)
library(dplyr)
library(RColorBrewer)

make.filename <- function(path,simnum,pi,mu,R,n_microbes,type)
{
  return(sprintf("%sR=%s_mu=%s_n_microbes=%s_out=%s_pi=%s_simnum=%s.csv",
                 path,R,mu,n_microbes,type,pi,simnum))
}

#pis and mus should be vectors of all the pi or mu values for this plot
create_edge_plots <- function(simnum,pis,mus,R,n_microbes,inpath,outpath)
{
  #flnm <- make.filename(inpath,simnum,pis[1],mus[2],R,n_microbes,"edges")
  edges <- data.frame()
  
  for(i.pi in pis)
  {
    for(j.mu in mus)
    {
      flnm <- make.filename(inpath,simnum,i.pi,j.mu,R,n_microbes,"edges")
      print(flnm)
      edges <- rbind(edges, read.csv(flnm))
    }
  }
  
  plt1 <- ggplot(edges, aes(x_microbe,y_microbe,fill=abs(mean))) + geom_tile() + 
    scale_y_reverse(breaks=1:29, labels=1:29, minor_breaks = seq(0.5,29.5),
                    position="right") + 
    scale_x_continuous(breaks=2:30, labels=2:30, position="top", 
                       minor_breaks = seq(1.5,30.5)) +
    scale_fill_gradient(low="#FFFFFF", high="#FF0000") +
    labs(fill="Coef.") + xlab("") + ylab("") +
    theme(panel.background=element_rect(fill="#FFFFFF"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_line(size=0.15,linetype='solid',
                                          color='black'))+
    facet_grid(pi ~ mu, switch="both",
               labeller = labeller(
                 mu = c(`0.8` = "mu=0.8", `1.6` = "mu=1.6"),
                 pi = c(`0.3` = "pi=0.3", `0.8` = "pi=0.8")
               ))
  ggsave(sprintf("%sR=%s_posterior_edges.png",outpath,R),plot=plt1,width=11,height=9.5)
  
  plt2 <- ggplot(edges, aes(x_microbe,y_microbe,fill=abs(true_B))) + geom_tile() + 
    scale_y_reverse(breaks=1:29, labels=1:29, minor_breaks = seq(0.5,29.5),
                    position="right") + 
    scale_x_continuous(breaks=2:30, labels=2:30, position="top", 
                       minor_breaks = seq(1.5,30.5)) +
    scale_fill_gradient(low="#FFFFFF", high="#FF0000") +
    labs(fill="Coef.") + xlab("") + ylab("") +
    theme(panel.background=element_rect(fill="#FFFFFF"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_line(size=0.15,linetype='solid',
                                          color='black'))+
    facet_grid(pi ~ mu, switch="both",
               labeller = labeller(
                 mu = c(`0.8` = "mu=0.8", `1.6` = "mu=1.6"),
                 pi = c(`0.3` = "pi=0.3", `0.8` = "pi=0.8")
               ))
  
  ggsave(sprintf("%sR=%s_true_edges.png",outpath,R),plot=plt2,width=11,height=9.5)
} 

create_MSE_plots <- function(simnum,pis,mus,R,n_microbes,inpath,outpath)
{
  mse <- data.frame()
  for(i.pi in pis)
  {
    for(j.mu in mus)
    {
      for(k.n in n_microbes)
      {
        flnm <- make.filename(inpath,simnum,i.pi,j.mu,R,k.n,"MSE")
        print(flnm)
        mse <- rbind(mse,read.csv(flnm))
      }
    }
  }
  plt1 <- ggplot(mse, aes(x=n_microbes, y=MSE, group=interaction(pi,mu),color=factor(pi), 
                  linetype=factor(mu)))  + 
    geom_line() + labs(x="Number of Microbes",color="pi value",linetype="mu value") +
    theme_bw() + scale_x_continuous(breaks=c(8,15,22)) + 
    scale_y_continuous(limits=c(0,ceiling(max(mse$MSE))),breaks=0:(ceiling(max(mse$MSE))*10)/10)
  ggsave(sprintf("%sR=%s_mse.png",outpath,R),plot=plt1)
}

create_node_plots <- function(simnum,pis,mus,R,n_microbes,inpath,outpath)
{
  nodes <- data.frame()
  for(i.pi in pis)
  {
    for(j.mu in mus)
    {
      flnm <- make.filename(inpath,simnum,pis[1],mus[2],R,n_microbes,"nodes")
      nodes <- rbind(nodes, read.csv(flnm))
    }
  }
}

create_interval_plots <- function()


create_plots <- function(simnum,pis,mus,R,n_microbes,inpath,outpath)
{
  for(r in R)
  {
    for(n_m in n_microbes)
    {
      create_edge_plots(simnum,pis,mus,r,n_m,inpath,outpath)
    }
    create_MSE_plots(simnum,pis,mus,r,n_microbes,inpath,outpath)
  }
}
create_plots(1,c(0.3,0.8),c(0.8,1.6),c(5,9),c(8,15,22),"results/simulation/","plots/simulation/")


