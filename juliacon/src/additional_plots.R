library(ggplot2)
library(dplyr)

# read in results
xis.1 <- read.csv("juliacon/results/simulation1_case1.csv")
xis.2 <- read.csv("juliacon/results/simulation1_case2.csv")
xis.3 <- read.csv("juliacon/results/simulation1_case3.csv")
xis.4 <- read.csv("juliacon/results/simulation1_case4.csv")
xis.5 <- read.csv("juliacon/results/simulation1_case5.csv")
xis.6 <- read.csv("juliacon/results/simulation1_case6.csv")
xis.7 <- read.csv("juliacon/results/simulation1_case7.csv")
xis.8 <- read.csv("juliacon/results/simulation1_case8.csv")

xis.1$microbe <- 1:30
xis.2$microbe <- 1:30
xis.3$microbe <- 1:30
xis.4$microbe <- 1:30
xis.5$microbe <- 1:30
xis.6$microbe <- 1:30
xis.7$microbe <- 1:30
xis.8$microbe <- 1:30

list.xis <- list(xis.1,xis.2,xis.3,xis.4,xis.5,xis.6,xis.7,xis.8)

case <- 0
for(xi in list.xis)
{
  case <- case + 1
  plt <- ggplot(data=xi,aes(x=microbe,y=Xi.posterior,fill=TrueXi)) + 
    geom_bar(stat="Identity") + xlab("Microbe") + ylab("Probability of influence") + 
    theme_minimal() + scale_fill_manual(values=c("#FF0000","#0000FF")) + 
    #labs(fill="Influential \n Microbe") + 
    guides(fill=FALSE) +
    scale_x_continuous(breaks=1:30,labels=1:30) + 
    ggtitle(sprintf("Influential Microbes: Simulation 1 Case %d",case)) +
    theme(plot.title=element_text(hjust=0.5))
  
  ggsave(sprintf("juliacon/plots/simulation1_case%d_posterior_xi.png",case),plot=plt)
}

### I'm not gonna do this the way I should because I don't want to think about how
### to loop through data frames


plt <- ggplot(data=xis.2,aes(x=microbe,y=Xi.posterior,fill=TrueXi)) + 
  geom_bar(stat="Identity") + xlab("Microbe") + ylab("Probability of influence") + 
  theme_minimal() + scale_fill_manual(values=c("#FF0000","#0000FF")) + 
  labs(fill="Influential 
 Microbe") + scale_x_continuous(breaks=1:30,labels=1:30)

ggsave("juliacon/plots/simulation1_case1",plot=plt)

plt <- ggplot(data=xis.3,aes(x=microbe,y=Xi.posterior,fill=TrueXi)) + 
  geom_bar(stat="Identity") + xlab("Microbe") + ylab("Probability of influence") + 
  theme_minimal() + scale_fill_manual(values=c("#FF0000","#0000FF")) + 
  labs(fill="Influential 
 Microbe") + scale_x_continuous(breaks=1:30,labels=1:30)

ggsave("juliacon/plots/simulation1_case1",plot=plt)

plt <- ggplot(data=xis.4,aes(x=microbe,y=Xi.posterior,fill=TrueXi)) + 
  geom_bar(stat="Identity") + xlab("Microbe") + ylab("Probability of influence") + 
  theme_minimal() + scale_fill_manual(values=c("#FF0000","#0000FF")) + 
  labs(fill="Influential 
 Microbe") + scale_x_continuous(breaks=1:30,labels=1:30)

ggsave("juliacon/plots/simulation1_case1",plot=plt)

plt <- ggplot(data=xis.5,aes(x=microbe,y=Xi.posterior,fill=TrueXi)) + 
  geom_bar(stat="Identity") + xlab("Microbe") + ylab("Probability of influence") + 
  theme_minimal() + scale_fill_manual(values=c("#FF0000","#0000FF")) + 
  labs(fill="Influential 
 Microbe") + scale_x_continuous(breaks=1:30,labels=1:30)

ggsave("juliacon/plots/simulation1_case1",plot=plt)

plt <- ggplot(data=xis.6,aes(x=microbe,y=Xi.posterior,fill=TrueXi)) + 
  geom_bar(stat="Identity") + xlab("Microbe") + ylab("Probability of influence") + 
  theme_minimal() + scale_fill_manual(values=c("#FF0000","#0000FF")) + 
  labs(fill="Influential 
 Microbe") + scale_x_continuous(breaks=1:30,labels=1:30)

ggsave("juliacon/plots/simulation1_case1",plot=plt)

plt <- ggplot(data=xis.7,aes(x=microbe,y=Xi.posterior,fill=TrueXi)) + 
  geom_bar(stat="Identity") + xlab("Microbe") + ylab("Probability of influence") + 
  theme_minimal() + scale_fill_manual(values=c("#FF0000","#0000FF")) + 
  labs(fill="Influential 
 Microbe") + scale_x_continuous(breaks=1:30,labels=1:30)

ggsave("juliacon/plots/simulation1_case1",plot=plt)

plt <- ggplot(data=xis.8,aes(x=microbe,y=Xi.posterior,fill=TrueXi)) + 
  geom_bar(stat="Identity") + xlab("Microbe") + ylab("Probability of influence") + 
  theme_minimal() + scale_fill_manual(values=c("#FF0000","#0000FF")) + 
  labs(fill="Influential 
 Microbe") + scale_x_continuous(breaks=1:30,labels=1:30)

ggsave("juliacon/plots/simulation1_case1",plot=plt)