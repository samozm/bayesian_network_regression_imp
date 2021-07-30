
library(ggplot2)
library(dplyr)
library(RColorBrewer)

make.filename <- function(path,simnum,pi,mu,R,nu,n_microbes,type,simtype=NULL)
{
  if(is.null(simtype))
  {
    return(sprintf("%sR=%s_mu=%s_n_microbes=%s_nu=%s_out=%s_pi=%s_simnum=%s.csv",
                   path,R,mu,n_microbes,nu,type,pi,simnum))
  }
  return(sprintf("%sR=%s_mu=%s_n_microbes=%s_nu=%s_out=%s_pi=%s_simnum=%s_type=%s.csv",
                path,R,mu,n_microbes,nu,type,pi,simnum,simtype))
  
}

plot_CIs <- function(edges)
{
  plt <- edges %>% filter(edge == 1) %>% ggplot() + geom_line(aes(x=`95CI`,
                                                                  y=edge))
  for(i in edges$edge)
  {
    plt <- plt + geom_line(aes(data=filter(edges,edge==i),x=`95CI`,y=edge))
  }
  return(plt)
}

#pis and mus should be vectors of all the pi or mu values for this plot
create_edge_plots <- function(simnum,pis,mus,R,nu,n_microbes,inpath,outpath,simtype=NULL)
{
  edges <- data.frame()
  
  for(i.pi in pis)
  {
    for(j.mu in mus)
    {
      flnm <- make.filename(inpath,simnum,i.pi,j.mu,R,nu,n_microbes,"edges",simtype)
      edges <- rbind(edges, read.csv(flnm))
    }
  }
  
  n <- length(edges$mean) / (length(pis)*length(mus))
  
  plt1 <- ggplot(edges, aes(x_microbe,y_microbe,fill=abs(mean))) + geom_tile() + 
    scale_y_reverse(breaks=1:(n-1), labels=1:(n-1), minor_breaks = seq(0.5,n-0.5),
                    position="right") + 
    scale_x_continuous(breaks=2:n, labels=2:n, position="top", 
                       minor_breaks = seq(1.5,n + 0.5)) +
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
               ),
               as.table = FALSE)
  if(is.null(simtype))
  {
    out_fl <- sprintf("%s%sR=%s_n_microbes=%s_nu=%s_posterior_edges.png",outpath,"edges/",
                      R,n_microbes,nu)
  }
  else {
    out_fl <- sprintf("%s%sR=%s_n_microbes=%s_nu=%s_type=%s_posterior_edges.png",
                      outpath,"edges/",R,n_microbes,nu,simtype)
  }
  
  ggsave(out_fl,plot=plt1,width=11,height=9.5)
  
  plt2 <- ggplot(edges, aes(x_microbe,y_microbe,fill=abs(true_B))) + geom_tile() + 
    scale_y_reverse(breaks=1:(n-1), labels=1:(n-1), minor_breaks = seq(0.5,n-0.5),
                    position="right") + 
    scale_x_continuous(breaks=2:n, labels=2:n, position="top", 
                       minor_breaks = seq(1.5,n+0.5)) +
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
               ),
               as.table = FALSE)
  if(is.null(simtype))
  {
    out_fl <- sprintf("%s%sR=%s_n_microbes=%s_nu=%s_true_edges.png",outpath,"edges/",
                      R,n_microbes,nu)
  } else {
    out_fl <- sprintf("%s%sR=%s_n_microbes=%s_nu=%s_type=%s_true_edges.png",outpath,"edges/",
                      R,n_microbes,nu,simtype)
  }
  ggsave(out_fl,plot=plt2,width=11,height=9.5)
}

create_MSE_plots <- function(simnum,pis,mus,R,nu,n_microbes,inpath,outpath,simtype=NULL)
{
  mse <- data.frame()
  for(i.pi in pis)
  {
    for(j.mu in mus)
    {
      for(k.n in n_microbes)
      {
        flnm <- make.filename(inpath,simnum,i.pi,j.mu,R,nu,k.n,"MSE",simtype)
        mse <- rbind(mse,read.csv(flnm))
      }
    }
  }
  plt1 <- ggplot(mse, aes(x=n_microbes, y=MSE, group=interaction(pi,mu),color=factor(pi), 
                  linetype=factor(mu)))  + 
    geom_line() + labs(x="Number of Microbes",color="pi value",linetype="mu value") +
    theme_bw() + scale_x_continuous(breaks=n_microbes) + 
    scale_y_continuous(limits=c(0,ceiling(max(mse$MSE))),
                       breaks=seq(0,ceiling(max(mse$MSE)),
                                  ceiling(max(mse$MSE))/10))
  
  if(is.null(simtype))
  {
    out_fl <- sprintf("%s%sR=%s_nu=%s_mse.png",outpath,"mse/",R,nu)
  } else {
    out_fl <- sprintf("%s%sR=%s_nu=%s_type=%s_mse.png",outpath,"mse/",R,nu,simtype)
  }
  ggsave(out_fl,plot=plt1)
}

create_node_plots <- function(simnum,pis,mus,R,nu,n_microbes,inpath,outpath,simtype=NULL)
{
  nodes <- data.frame()
  for(i.pi in pis)
  {
    for(j.mu in mus)
    {
      flnm <- make.filename(inpath,simnum,i.pi,j.mu,R,nu,n_microbes,"nodes",simtype)
      nodes <- rbind(nodes, read.csv(flnm))
    }
  }
  
  n <- length(nodes$TrueXi) / (length(pis)*length(mus))
  
  nodes$microbe <- rep(1:n,length(pis)*length(mus))
  plt <- ggplot(data=nodes,aes(x=microbe,y=Xi.posterior,fill=TrueXi)) + 
    geom_bar(stat="Identity") + xlab("Microbe") + ylab("Probability of influence") + 
    theme_minimal() +
    scale_fill_manual(values=c("#FF0000","#0000FF")) + 
    guides(fill="none") +
    scale_x_continuous(breaks=1:n,labels=1:n,position="top") + 
    theme(plot.title=element_text(hjust=0.5), panel.background = element_rect(fill="#FFFFFF")) +
    scale_y_continuous(position="right") + 
    facet_grid(pi ~ mu, switch="both",
              labeller = labeller(
                mu = c(`0.8` = "mu=0.8", `1.6` = "mu=1.6"),
                pi = c(`0.3` = "pi=0.3", `0.8` = "pi=0.8")
              ),
              as.table = FALSE)
  if(is.null(simtype))
  {
    out_fl <- sprintf("%s%sR=%s_n_microbes=%s_nu=%s_nodes.png",outpath,"nodes/",R,
                      n_microbes,nu)
  } else {
    out_fl <- sprintf("%s%sR=%s_n_microbes=%s_nu=%s_type=%s_nodes.png",outpath,"nodes/",R,
                      n_microbes,nu,simtype)
  }
  ggsave(out_fl,plot=plt,width=11,height=9.5)
}

create_interval_plots <- function(simnum,pi,mu,R,nu,n_microbes,inpath,outpath,simtype=NULL)
{
  flnm <- make.filename(inpath,simnum,pi,mu,R,nu,n_microbes,"edges",simtype)
  edges <- read.csv(flnm)
  
  n <- length(edges$true_B)
  
  edges$edge <- 1:n
  edges <- transform(edges,rej=ifelse(X0.025 > 0 | X0.975 < 0,TRUE,FALSE))
  edges <- transform(edges,nonzero_mean=ifelse(true_B != 0.0,TRUE,FALSE))
  plt1 <- edges %>% ggplot() + geom_errorbar(aes(x=factor(edge),ymin=X0.025,
                                                                 ymax=X0.975,
                                                                 color=rej)) +
          theme_bw() + theme(panel.grid.major = element_blank()) + labs(color="Edge Appears Influential") +
          #scale_x_continuous(limits=c(1,n)) + coord_flip() +
          scale_y_continuous(limits=c(-15,15)) +
          coord_flip() + scale_x_discrete(labels=NULL,expand=expansion(add=4)) +
          geom_hline(aes(yintercept=mu[1]),linetype="dashed",color="blue") + 
          geom_point(aes(x=factor(edge),y=mean),shape=5) + ylab("Edge Effect") +
          geom_hline(aes(yintercept=0),linetype="dashed",color="green") +
          facet_grid(nonzero_mean ~.,scales="free_y",space="free_y", labeller = labeller(
            nonzero_mean = c(`TRUE` = "True mean nonzero", `FALSE` = "True mean zero")
          )) 
  if(is.null(simtype)){
    out_fl <- sprintf("%s%sR=%s_pi=%s_mu=%s_n_microbes=%s_nu=%s_intervals.png",
                      outpath,"intervals/",R,pi,mu,n_microbes,nu)
  } else {
    out_fl <- sprintf("%s%sR=%s_pi=%s_mu=%s_n_microbes=%s_nu=%s_type=%s_intervals.png",
                      outpath,"intervals/",R,pi,mu,n_microbes,nu,simtype)
  }
  ggsave(out_fl,plot=plt1,width=11,height=18)
  
}
#create_interval_plots(1,0.3,0.8,5,10,8,"results/simulation/","plots/simulation/")


create_plots <- function(simnum,pis,mus,R,nus,n_microbes,inpath,outpath,simtypes=c(NULL))
{
  for(s in simtypes)
  { 
    for(r.i in 1:length(R))
    {
      for(n_m in n_microbes)
      {
        create_edge_plots(simnum,pis,mus,R[r.i],nus[r.i],n_m,inpath,outpath,s)
        create_node_plots(simnum,pis,mus,R[r.i],nus[r.i],n_m,inpath,outpath,s)
        for(pi in pis)
        {
          for(mu in mus)
          {
            create_interval_plots(simnum,pi,mu,R[r.i],nus[r.i],n_m,inpath,outpath,s)
          }
        }
      }
      create_MSE_plots(simnum,pis,mus,R[r.i],nus[r.i],n_microbes,inpath,outpath,s)
    }
  }
}
create_plots(1,c(0.3,0.8),c(0.8,1.6),c(5,10,15),c(10,15,20),c(8,15,22),"results/simulation/unrealistic/","plots/simulation/unrealistic/")
create_plots(1,c(0.3,0.8),c(0.8,1.6),c(15,15,25),c(16,17,30),c(8,15,22),"results/simulation/unrealistic/","plots/simulation/unrealistic/")
create_plots(1,c(0.3,0.8),c(0.8,1.6),c(5),c(7),c(8,15,22),"results/simulation/unrealistic/","plots/simulation/unrealistic/")
create_plots(1,c(0.3,0.8),c(0.8,1.6),c(9),c(10),c(8,15,22),"results/simulation/unrealistic/","plots/simulation/unrealistic/")

create_plots(2,c(0.3,0.8),c(0.8,1.6),c(9),c(10),c(8,22),"results/simulation/realistic/",
             "plots/simulation/realistic/",
             c("additive_phylo", "additive_random", "interaction_phylo", 
               "interaction_random", "redundant_phylo", "redundant_random"))

