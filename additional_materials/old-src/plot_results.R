
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(stringr)

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
create_edge_plots <- function(simnum,pis,mus,R,nu,n_microbes,inpath,outpath,samplesize,simtype=NULL,pref="",edge_mu=NULL)
{
  edges <- data.frame()
  
  for(i.pi in pis)
  {
    for(j.mu in mus)
    {
      flnm <- make.filename(inpath,simnum,i.pi,j.mu,R,nu,n_microbes,"edges",samplesize,simtype,edge_mu)
      print(flnm)
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
                                          color='black'),
          plot.title = element_text(hjust = 0.5))+
    ggtitle("Posterior Edge Coefficients")+
    facet_grid(pi ~ mu, switch="both",
               labeller = labeller(
                 mu = c(`0.8` = "mu=0.8", `1.6` = "mu=1.6"),
                 pi = c(`0` = "pi=0.0", `0.3` = "pi=0.3", `0.8` = "pi=0.8")
               ),
               as.table = FALSE)
  if(is.null(simtype))
  {
    out_fl <- sprintf("%s%s%sR=%s_n_microbes=%s_nu=%s_samplesize=%s_posterior_edges.png",outpath,
                      "edges/",pref,
                      R,n_microbes,nu,samplesize)
  }
  else {
    out_fl <- sprintf("%s%s%sR=%s_edge_mu=%s_n_microbes=%s_nu=%s_samplesize=%s_type=%s_posterior_edges.png",
                      outpath,"edges/",pref,R,edge_mu,
                      n_microbes,nu,samplesize,simtype)
  }
  print(out_fl)
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
                                          color='black'),
          plot.title = element_text(hjust = 0.5))+
    ggtitle("True Edge Coefficients")+
    facet_grid(pi ~ mu, switch="both",
               labeller = labeller(
                 mu = c(`0.8` = "mu=0.8", `1.6` = "mu=1.6"),
                 pi = c(`0` = "pi=0.0", `0.3` = "pi=0.3", `0.8` = "pi=0.8")
               ),
               as.table = FALSE)
  if(is.null(simtype))
  {
    out_fl <- sprintf("%s%s%sR=%s_n_microbes=%s_nu=%s_samplesize=%s_true_edges.png",outpath,
                      "edges/",pref,
                      R,n_microbes,nu,samplesize)
  } else {
    out_fl <- sprintf("%s%s%sR=%s_edge_mu=%s_n_microbes=%s_nu=%s_samplesize=%s_type=%s_true_edges.png",
                      outpath,"edges/",pref,
                      R,edge_mu,n_microbes,
                      nu,samplesize,simtype)
  }
  ggsave(out_fl,plot=plt2,width=11,height=9.5)
}

create_MSE_plots <- function(simnum,pis,mus,R,nu,n_microbes,inpath,outpath,samplesize,simtype=NULL,pref="",edge_mu=NULL)
{
  mse <- data.frame()
  for(i.pi in pis)
  {
    for(j.mu in mus)
    {
      for(k.n in n_microbes)
      {
        flnm <- make.filename(inpath,simnum,i.pi,j.mu,R,nu,k.n,"MSE",samplesize,
                              simtype,edge_mu)
        mse <- rbind(mse,read.csv(flnm))
      }
    }
  }
  plt1 <- ggplot(mse, aes(x=n_microbes, y=MSE, group=interaction(pi,mu),color=factor(pi), 
                  linetype=factor(mu)))  + 
    geom_line() + labs(x="Number of Microbes",color="pi value",linetype="mu value",
                       title="MSE - Coefficients") +
    theme_bw() + scale_x_continuous(breaks=n_microbes) + 
    scale_y_continuous(limits=c(0,ceiling(max(mse$MSE))),
                       breaks=seq(0,ceiling(max(mse$MSE)),
                                  ceiling(max(mse$MSE))/10))
  
  if(is.null(simtype))
  {
    out_fl <- sprintf("%s%s%sR=%s_nu=%s_samplesize=%s_mse.png",outpath,"mse/",
                      pref,R,nu,samplesize)
  } else {
    out_fl <- sprintf("%s%s%sR=%s_edge_mu=%s_nu=%s_samplesize=%s_type=%s_mse.png",outpath,
                      "mse/",pref,R,edge_mu,nu,samplesize,simtype)
  }
  ggsave(out_fl,plot=plt1)
  
  plt2 <- ggplot(mse, aes(x=n_microbes, y=MSEy, group=interaction(pi,mu),color=factor(pi), 
                          linetype=factor(mu)))  + 
    geom_line() + labs(x="Number of Microbes",color="pi value",linetype="mu value",
                       title="MSE - Response") +
    theme_bw() + scale_x_continuous(breaks=n_microbes) + 
    scale_y_continuous(limits=c(0,ceiling(max(mse$MSEy))),
                       breaks=seq(0,ceiling(max(mse$MSEy)),
                                  ceiling(max(mse$MSEy))/10))
  
  if(is.null(simtype))
  {
    out_fl <- sprintf("%s%s%sR=%s_nu=%s_samplesize=%s_mse-y.png",outpath,"mse/",
                      pref,R,nu,samplesize)
  } else {
    out_fl <- sprintf("%s%s%sR=%s_edge_mu=%s_nu=%s_samplesize=%s_type=%s_mse-y.png",outpath,
                      "mse/",pref,R,edge_mu,nu,samplesize,simtype)
  }
  ggsave(out_fl,plot=plt2)
}

create_node_plots <- function(simnum,pis,mus,R,nu,n_microbes,inpath,outpath,samplesize,simtype=NULL,pref="",edge_mu=NULL)
{
  nodes <- data.frame()
  for(i.pi in pis)
  {
    for(j.mu in mus)
    {
      flnm <- make.filename(inpath,simnum,i.pi,j.mu,R,nu,n_microbes,"nodes",
                            samplesize,simtype,edge_mu)
      print(flnm)
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
    scale_y_continuous(position="right",breaks = c(0,0.5,1),limits=c(0,1)) + 
    facet_grid(pi ~ mu, switch="both",
              labeller = labeller(
                mu = c(`0.8` = "mu=0.8", `1.2` = "mu=1.2", `1.4` = "mu=1.4", `1.6` = "mu=1.6"),
                pi = c(`0` = "pi=0.0", `0.3` = "pi=0.3", `0.5` = "pi=0.5", `0.7` = "pi=0.7", `0.8` = "pi=0.8")
              ),
              as.table = FALSE)
  if(is.null(simtype))
  {
    out_fl <- sprintf("%s%s%sR=%s_n_microbes=%s_nu=%s_samplesize=%s_nodes.png",outpath,
                      "nodes/",pref,R,
                      n_microbes,nu,samplesize)
  } else {
    out_fl <- sprintf("%s%s%sR=%s_edge_mu=%s_n_microbes=%s_nu=%s_samplesize=%s_type=%s_nodes.png",
                      outpath,"nodes/",pref,R,edge_mu,
                      n_microbes,nu,samplesize,simtype)
  }
  ggsave(out_fl,plot=plt,width=11,height=9.5)
}

create_interval_plots <- function(simnum,pi,mu,R,nu,n_microbes,inpath,outpath,samplesize,simtype=NULL,pref="",edge_mu=NULL)
{
  flnm <- make.filename(inpath,simnum,pi,mu,R,nu,n_microbes,"edges",samplesize,simtype,edge_mu)
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
          geom_point(aes(x=factor(edge),y=mean),shape=5) + ylab("Edge Effect") +
          geom_hline(aes(yintercept=0),linetype="dashed",color="green") +
          facet_grid(nonzero_mean ~.,scales="free_y",space="free_y", labeller = labeller(
            nonzero_mean = c(`TRUE` = "True mean nonzero", `FALSE` = "True mean zero")
          )) 
  
  plt2 <- ggplot(edges, aes(x_microbe,y_microbe,fill=factor(rej),width=0.98,height=0.98)) + geom_tile() + 
    scale_y_reverse(breaks=1:(n-1), labels=1:(n-1), minor_breaks = seq(0.5,n-0.5),
                    position="right") + 
    scale_x_continuous(breaks=2:n, labels=2:n, position="top", 
                       minor_breaks = seq(1.5,n + 0.5)) +
    scale_fill_discrete(type=c("#FFFFFF","#FF0000"),labels=c("","Influential")) +
    labs(fill="Influential Edge") + xlab("") + ylab("") +
    theme(panel.background=element_rect(fill="#FFFFFF"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_line(size=0.15,linetype='solid',
                                          color='black'),
          plot.title = element_text(hjust = 0.5))+
    ggtitle("Posterior Edge Coefficients")+
    facet_grid(pi ~ mu, switch="both",
               labeller = labeller(
                 mu = c(`0.8` = "mu=0.8", `1.6` = "mu=1.6"),
                 pi = c(`0` = "pi=0.0", `0.3` = "pi=0.3", `0.8` = "pi=0.8")
               ),
               as.table = FALSE)
  
  if(is.null(simtype)){
    out_fl <- sprintf("%s%s%sR=%s_pi=%s_mu=%s_n_microbes=%s_nu=%s_samplesize=%s_intervals.png",
                      outpath,"intervals/",pref,R,pi,mu,n_microbes,nu,samplesize)
    plt1 <- plt1 + geom_hline(aes(yintercept=mu[1]),linetype="dashed",color="blue")
    out_fl2 <- sprintf("%s%s%sR=%s_pi=%s_mu=%s_n_microbes=%s_nu=%s_samplesize=%s_edge_grid.png",
                      outpath,"intervals/",pref,R,pi,mu,n_microbes,nu,samplesize)
  } else {
    out_fl <- sprintf("%s%s%sR=%s_edge_mu=%s_pi=%s_mu=%s_n_microbes=%s_nu=%s_samplesize=%s_type=%s_intervals.png",
                      outpath,"intervals/",pref,R,edge_mu,pi,mu,n_microbes,nu,samplesize,simtype)
    plt1 <- plt1 + geom_hline(aes(yintercept=edge_mu),linetype="dashed",color="blue")
    out_fl2 <- sprintf("%s%s%sR=%s_edge_mu=%s_pi=%s_mu=%s_n_microbes=%s_nu=%s_samplesize=%s_type=%s_edge_grid.png",
                      outpath,"intervals/",pref,R,edge_mu,pi,mu,n_microbes,nu,samplesize,simtype)
  }
  ggsave(out_fl,plot=plt1,width=11,height=18)
  ggsave(out_fl2,plot=plt2,width=11,height=9.5)
  
}

create_fr_plots <- function(simnum,pis,mus,R,nu,n_microbes,inpath,outpath,samplesize,simtype=NULL,pref="",edge_mu=NULL)
{
  
  n <- length(n_microbes) * length(pis) * length(mus)
  fr.e <- data.frame(n_microbes=n_microbes,fpr_edges=rep(0,n), fpr_nodes=rep(0,n),
                      fnr_edges=rep(0,n), fnr_nodes=rep(0,n),
                      pi=rep(0,n),mu=rep(0,n))
  l <- 1
  
  for(i.pi in pis)
  {
    for(j.mu in mus)
    {
      for(k.n in n_microbes)
      {
        
        flnm <- make.filename(inpath,simnum,i.pi,j.mu,R,nu,k.n,"edges",samplesize,
                              simtype,edge_mu)
        edges <- read.csv(flnm)
        
        o <- length(edges$true_B)
        edges$edge <- 1:o
        edges <- transform(edges,rej=ifelse(X0.025 > 0 | X0.975 < 0,1,0))
        edges <- transform(edges,nonzero_mean=ifelse(true_B != 0.0,1,0))
        edges <- transform(edges,zero_mean=ifelse(true_B == 0.0,1,0))
        edges <- transform(edges,fp=ifelse(rej==1 & nonzero_mean==0,1,0))
        edges <- transform(edges,fn=ifelse(rej==0 & nonzero_mean==1,1,0))
        
        fr.e$fpr_edges[l] <- sum(edges$fp) / sum(1 - edges$nonzero_mean)
        fr.e$fnr_edges[l] <- sum(edges$fn) / sum(edges$nonzero_mean)
        fr.e$pi[l] <- i.pi
        fr.e$mu[l] <- j.mu
        fr.e$n_microbes[l] <- k.n
        
        flnm2 <- make.filename(inpath,simnum,i.pi,j.mu,R,nu,k.n,"nodes",samplesize,
                               simtype,edge_mu)
        nodes <- read.csv(flnm2)
        
        m <- length(nodes$TrueXi)
        nodes$node <- 1:m
        nodes <- transform(nodes, rej=ifelse(Xi.posterior>=0.5,1,0))
        nodes <- transform(nodes, TrueXi=ifelse(TrueXi=='true',TRUE,FALSE))
        nodes <- transform(nodes,fp=ifelse(rej==1 & TrueXi==FALSE,1,0))
        nodes <- transform(nodes,fn=ifelse(rej==0 & TrueXi==TRUE,1,0))
        
        fr.e$fpr_nodes[l] <- sum(nodes$fp) / sum(1 - nodes$TrueXi)
        fr.e$fnr_nodes[l] <- sum(nodes$fn) / sum(nodes$TrueXi)
        
        l <- l+1
      }
    }
  }
  
  plt1 <- ggplot(fr.e, aes(x=n_microbes, y=fpr_edges, group=interaction(pi,mu),color=factor(pi), 
                          linetype=factor(mu)))  + 
  geom_line() + labs(x="Number of Microbes",color="pi value",linetype="mu value",
                     title="False Positive Rate - Edges", y = "False Positive Rate") +
  theme_bw() + scale_x_continuous(breaks=n_microbes) + 
  scale_y_continuous(limits=c(0,1),
                     breaks=seq(0,1,
                                1/10))
  plt2 <- ggplot(fr.e, aes(x=n_microbes, y=fpr_nodes, group=interaction(pi,mu),color=factor(pi), 
                            linetype=factor(mu)))  + 
    geom_line() + labs(x="Number of Microbes",color="pi value",linetype="mu value",
                       title="False Positive Rate - Nodes",  y = "False Positive Rate") +
    theme_bw() + scale_x_continuous(breaks=n_microbes) + 
    scale_y_continuous(limits=c(0,1),
                       breaks=seq(0,1,
                                  1/10))
  
  plt3 <- ggplot(fr.e, aes(x=n_microbes, y=fnr_edges, group=interaction(pi,mu),color=factor(pi), 
                            linetype=factor(mu)))  + 
    geom_line() + labs(x="Number of Microbes",color="pi value",linetype="mu value",
                       title="False Negative Rate - Edges", y = "False Negative Rate") +
    theme_bw() + scale_x_continuous(breaks=n_microbes) + 
    scale_y_continuous(limits=c(0,1),
                       breaks=seq(0,1,
                                  1/10))
  plt4 <- ggplot(fr.e, aes(x=n_microbes, y=fnr_nodes, group=interaction(pi,mu),color=factor(pi), 
                            linetype=factor(mu)))  + 
    geom_line() + labs(x="Number of Microbes",color="pi value",linetype="mu value",
                       title="False Negative Rate - Nodes",  y = "False Negative Rate") +
    theme_bw() + scale_x_continuous(breaks=n_microbes) + 
    scale_y_continuous(limits=c(0,1),
                       breaks=seq(0,1,
                                  1/10))
  
  if(is.null(simtype))
  {
    out_fl1 <- sprintf("%s%s%sR=%s_nu=%s_samplesize=%s_edge_fpr.png",outpath,"fr/",
                      pref,R,nu,samplesize)
    out_fl2 <- sprintf("%s%s%sR=%s_nu=%s_samplesize=%s_node_fpr.png",outpath,"fr/",
                       pref,R,nu,samplesize)
    out_fl3 <- sprintf("%s%s%sR=%s_nu=%s_samplesize=%s_edge_fnr.png",outpath,"fr/",
                       pref,R,nu,samplesize)
    out_fl4 <- sprintf("%s%s%sR=%s_nu=%s_samplesize=%s_node_fnr.png",outpath,"fr/",
                       pref,R,nu,samplesize)
  } else {
    out_fl1 <- sprintf("%s%s%sR=%s_edge_mu=%s_nu=%s_samplesize=%s_type=%s_edge_fpr.png",outpath,
                      "fr/",pref,R,edge_mu,nu,samplesize,simtype)
    out_fl2 <- sprintf("%s%s%sR=%s_edge_mu=%s_nu=%s_samplesize=%s_type=%s_node_fpr.png",outpath,
                       "fr/",pref,R,edge_mu,nu,samplesize,simtype)
    out_fl3 <- sprintf("%s%s%sR=%s_edge_mu=%s_nu=%s_samplesize=%s_type=%s_edge_fnr.png",outpath,
                       "fr/",pref,R,edge_mu,nu,samplesize,simtype)
    out_fl4 <- sprintf("%s%s%sR=%s_edge_mu=%s_nu=%s_samplesize=%s_type=%s_node_fnr.png",outpath,
                       "fr/",pref,R,edge_mu,nu,samplesize,simtype)
  }
  ggsave(out_fl1,plot=plt1)
  ggsave(out_fl2,plot=plt2)
  ggsave(out_fl3,plot=plt3)
  ggsave(out_fl4,plot=plt4)
}
#create_fdr_plots(2,c(0.3,0.8),c(0.8,1.6),9,10,c(8,22),"results/simulation/realistic/",
#                 "plots/simulation/realistic/",100,
#                 "redundant_phylo","",0.4)


create_plots <- function(simnum,pis,mus,R,nus,n_microbes,inpath,outpath,samplesize,simtypes=c(NULL),pref="",edge_mu=NULL)
{
  if(!is.null(simtypes)){
    for(s in simtypes)
    { 
      plot_loops(simnum,pis,mus,R,nus,n_microbes,inpath,outpath,samplesize,s,pref,edge_mu)
    }
  }else{
    plot_loops(simnum,pis,mus,R,nus,n_microbes,inpath,outpath,samplesize,NULL,pref,edge_mu)
  }
}
plot_loops <- function(simnum,pis,mus,R,nus,n_microbes,inpath,outpath,samplesizes,simtype=NULL,pref="",edge_mu=NULL){
  for(samplesize in samplesizes)
    {
    for(r.i in 1:length(R))
    {
      for(n_m in n_microbes)
      {
        create_edge_plots(simnum,pis,mus,R[r.i],nus[r.i],n_m,inpath,outpath,samplesize,simtype,pref,edge_mu)
        create_node_plots(simnum,pis,mus,R[r.i],nus[r.i],n_m,inpath,outpath,samplesize,simtype,pref,edge_mu)
        for(pi in pis)
        {
          for(mu in mus)
          {
            create_interval_plots(simnum,pi,mu,R[r.i],nus[r.i],n_m,inpath,outpath,samplesize,simtype,pref,edge_mu)
          }
        }
      }
      create_MSE_plots(simnum,pis,mus,R[r.i],nus[r.i],n_microbes,inpath,outpath,samplesize,simtype,pref,edge_mu)
      create_fr_plots(simnum,pis,mus,R[r.i],nus[r.i],n_microbes,inpath,outpath,samplesize,simtype,pref,edge_mu)
    }
  }
}

# unrealistic, power, 100 samples
#create_plots(1,c(0.3,0.8),c(0.8,1.6),c(5,10,15),c(10,15,20),c(8,15,22),
#             "results/simulation/unrealistic/","plots/simulation/unrealistic/",
#             c(100))
#create_plots(1,c(0.3,0.8),c(0.8,1.6),c(15,15,25),c(16,17,30),c(8,15,22),
#             "results/simulation/unrealistic/","plots/simulation/unrealistic/",
#             c(100))
#create_plots(1,c(0.3,0.8),c(0.8,1.6),c(5),c(7),c(8,15,22),
#             "results/simulation/unrealistic/","plots/simulation/unrealistic/",
#             c(100))
#create_plots(1,c(0.3,0.8),c(0.8,1.6),c(9),c(10),c(8,15,22),
#             "results/simulation/unrealistic/","plots/simulation/unrealistic/",
#             c(100))

# unrealistic, power, 500 samples
#create_plots(1,c(0.3,0.8),c(0.8,1.6),c(9),c(10),c(8,15,22),
#             "results/simulation/unrealistic/","plots/simulation/unrealistic/",
#             c(500))

# CHTC unrealistic, null, 500 samples
# create_plots(1,c("0.0"),c(0.8,1.6),c(9),c(10),c(8,15,22),
#             "results/simulation/chtc-unrealistic/results/","plots/simulation/chtc-unrealistic/nullsim/",
#             c(500),NULL,"nullsim_")

## CHTC UNREALISTIC 500,100 samples
create_plots(1,c(0.3,0.8),c(0.8,1.6),c(7),c(10),c(8,15,22),
            "results/simulation/chtc-unrealistic-O0-R7/","results/plots/simulation/chtc-unrealistic-O0-R7/",
             c(500,100))

create_plots(1,c(0.3,0.8),c(0.8,1.6),c(9),c(10),c(8,15,22),
             "results/simulation/chtc-unrealistic-O0/","results/plots/simulation/chtc-unrealistic/",
             c(500,100))


## CHTC REALISTIC 500,1000 samples
#create_plots(2,c(0.3,0.8),c(0.8,1.6),c(9),c(10),c(8,22),
#             "results/simulation/chtc-realistic/results/","plots/simulation/chtc-realistic/",
#             c(500,1000),c("additive_phylo", "additive_random", "interaction_phylo", 
#                          "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)

create_plots(2,c(0.3,0.8),c(0.8,1.6),c(7),c(10),c(8,22),
             "results/simulation/chtc-realistic-O0-R7/","results/plots/simulation/chtc-realistic-O0-R7/",
             c(500,1000),c("additive_phylo", "additive_random", "interaction_phylo", 
                          "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)

#create_plots(1,c(0.8),c("2.0"),c(9),c(10),c(8),
#             "results/simulation/unrealistic/","plots/simulation/unrealistic/",
#             c(500))

# unrealistic, power, 1000 samples
#create_plots(1,c(0.3,0.8),c(0.8,1.6),c(9),c(10),c(8,15,22),
#             "results/simulation/unrealistic/","plots/simulation/unrealistic/",
#             c(1000))

#create_plots(2,c(0.3,0.8),c(0.8,1.6),c(9),c(10),c(8,22),"results/simulation/realistic/",
#             "plots/simulation/realistic/",c(100),
#             c("additive_phylo", "additive_random", "interaction_phylo", 
#               "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)


#create_plots(2,c(0.3,0.8),c(0.8,1.6),c(9),c(10),c(8,22),"results/simulation/realistic/",
#             "plots/simulation/realistic/",c(500),
#             c("additive_phylo", "additive_random", "interaction_phylo", 
#               "interaction_random", "redundant_phylo", "redundant_random"),"",0.4)


### intermediate mu and pi values 
#create_node_plots(1,c(0.3,0.5,0.7,0.8),c('1.0',1.2,1.4,1.6),9,10,22,
#                  "results/simulation/chtc-unrealistic/","plots/simulation/chtc-unrealistic/",
#                  c(500))

