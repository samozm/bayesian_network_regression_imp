
library(ggplot2,RColorBrewer)

#pis and mus should be vecotrs of all the pi or mu values for this plot
create_edge_plots <- function(simnum,pis,mus,R,n_microbes,inpath,outpath)
{
  flnm <- sprintf("%ssim%s_pi%s-mu%s-R%s-%smicrobes",
                  inpath,simnum,pi,mu,R,n_microbes)
  edges <- read.csv(sprintf("%s_edges.csv",flnm))
  nodes <- read.csv(sprintf("%s_nodes.csv",flnm))
  mse <- read.csv(sprintf("%s_MSE.csv",flnm))
  
  for(i.pi in c(0.1,0.8))
  {
    for(j.mu in c(0.8,1.6))
    {
      flnm <- sprintf("%ssim%s_pi%s-mu%s-R%s-%smicrobes",
                      inpath,simnum,i.pi,j.mu,R,n_microbes)
      edges <- rbind(edges,
                     read.csv(sprintf("%s_edges.csv",flnm)))
      nodes <- rbind(nodes,
                     read.csv(sprintf("%s_nodes.csv",flnm)))
      mse <- rbind(mse,
                   read.csv(sprintf("%s_MSE.csv",flnm)))
      
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
                 pi = c(`0.1` = "pi=0.1", `0.8` = "pi=0.8")
               ))
  ggsave(sprintf("%sposterior_edges.png",outpath),plot=plt1,width=11,height=9.5)
  
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
                 pi = c(`0.1` = "pi=0.1", `0.8` = "pi=0.8")
               ))
  
  ggsave(sprintf("%strue_edges.png",outpath),plot=plt2)
} 
create_plots <- function(simnum,pis,mus,R,n_microbes,inpath,outpath)
{
  create_edge_plots(simnum,pis,mus,R,n_microbes,inpath,outpath)
}
create_plots(1,c(0.1,0.8),c(0.8,1.6),5,8,"juliacon/results/","juliacon/plots/")

  
