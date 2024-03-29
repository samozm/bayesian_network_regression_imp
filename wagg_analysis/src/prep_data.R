library(data.table)
library(gdata)
library(readr)
library(ggplot2)

set.seed(122522)

setwd("~/Documents/masters/research/bayesian_network_regression/bayesian_network_regression_imp/wagg_analysis")

## RUN SPIEC-EASI
# ------------ --------------------------  ------------ #
# ------------ --------------------------  ------------ #
# Script for making a meta network of bacterial and fungi
# For details on the "spiec.easi" function see: 
# Kurtz, Z. D. et al. Sparse and Compositionally Robust Inference of Microbial Ecological Networks. PLoS Comput. Biol. 11 (2015), doi:10.1371/journal.pcbi.1004226
# To instal SpiecEasi package
library(devtools)
#install_github("zdk123/SpiecEasi")
# ------------ --------------------------  ------------ #
# ------------ --------------------------  ------------ #

# Load libraries
library(vegan)
library(SpiecEasi)
Func = read.csv("Wagg_etal_data/EcosystemFunctions.csv")
head(Func)
BACT  = read.csv("Wagg_etal_data/Bacterial_OTU_data.csv",row.names="X")
BACT = t(BACT)

# Checking the rows in the data and OTU matricies match...
rownames(BACT)[match(as.character(Func$ET), rownames(BACT))]
BACT = BACT[match(as.character(Func$ET), rownames(BACT)),]

# re-scale OTU counts to the proportion of the minimum sequencing depth
Bact = min(rowSums(BACT))*(BACT/rowSums(BACT))

# Remove OTUs that occur less than 40/45 times in the entire dataset
Bact40 = Bact[,which(specnumber(Bact, MARGIN=2) > 40)]

# Run the spiec.easi function (NOTE: this will take quie a while to run
se.est.40 = spiec.easi(Bact40, method='mb', lambda.min.ratio=1e-2, nlambda=50, icov.select.params=list(rep.num=100))

bm.40 <- symBeta(getOptBeta(se.est.40), mode="maxabs")
bm.40 <- abs(bm.40)

bm.mat.40 <- as.matrix(bm.40)
colnames(bm.mat.40) = rownames(bm.mat.40) = colnames(Bact40)

# Check matrix
bm.mat.40[1:20, 1:20]
write.csv(bm.mat.40,"Wagg_etal_data/Meta_network_weighted_40.csv") 


## 40 cutoff - 100 total samples

META_MATRIX <- read.csv("Wagg_etal_data/Meta_network_weighted_40.csv")

OTU_DATA <- read.csv("Wagg_etal_data/Bacterial_OTU_data.csv")
OTU_DATA[1:10,1:10]


OTU_DATAT <- transpose(OTU_DATA)
OTU_DATAT <- OTU_DATAT[-c(1),]

MATLEN <- dim(META_MATRIX)[1]

rownames(OTU_DATAT) <- colnames(OTU_DATA[,-c(1)])
colnames(OTU_DATAT) <- parse_number(OTU_DATA$X)

OTU_DATAT <- OTU_DATAT[,order(as.numeric(colnames(OTU_DATAT)))]
OTU_DATAT[1:10,1:10]

B.META <- META_MATRIX
dim(B.META)
B.META[1:10,1:10]

rownames(B.META) <- parse_number(gsub('\\.','',B.META$X))
B.META <- B.META[,c(-1)]
dim(B.META)
colnames(B.META) <- parse_number(gsub('\\.','',colnames(B.META)))
B.META[MATLEN-10:MATLEN,MATLEN-10:MATLEN]

B.META <- B.META[order(as.numeric(rownames(B.META))),order(as.numeric(colnames(B.META)))]
B.META[1:10,1:10]

upperTriangle(B.META,byrow=T)

colnames(Bact40) <- parse_number(gsub('\\.','',colnames(Bact40)))

OTU_DATAT_RES_TMP <- OTU_DATAT[,colnames(OTU_DATAT) %in% colnames(B.META)]
dim(OTU_DATAT_RES_TMP)

LEN_OTU <- dim(OTU_DATAT_RES_TMP)[1]
OTU_DATAT_RES <- data.frame(matrix(0,LEN_OTU*2,dim(OTU_DATAT_RES_TMP)[2]))
colnames(OTU_DATAT_RES) <- colnames(OTU_DATAT_RES_TMP)

func <- data.frame(matrix(0,LEN_OTU*2,dim(Func)[2]))
colnames(func) <- colnames(Func)

keep_prob <- 0.90

new_nleach_sd <- (sd(Func$Nleach[1:50]))^(1/50)
new_pleach_sd <- (sd(Func$Pleach[1:50]))/2
new_decom_sd <- (sd(Func$Decom[1:50]))^(1/50)

for(i in 1:LEN_OTU)
{
  num_necesary <- length(OTU_DATAT_RES_TMP[i,])
  OTU_DATAT_RES[i,] <- OTU_DATAT_RES_TMP[i,]
  for (j in 1:num_necesary)
  {
    OTU_DATAT_RES[i+LEN_OTU,j] <- ifelse(OTU_DATAT_RES_TMP[i,j] > 0, rbinom(1,1,keep_prob), rbinom(1,1,1-keep_prob))
  }
  func[i,] <- Func[i,]
  func[i+LEN_OTU,] <- Func[i,]
  
  new_nleach <- func[i+LEN_OTU,"Nleach"] + round(rnorm(1,0,new_nleach_sd),digits=2)
  if(new_nleach <= 0) { new_nleach <- rchisq(1,3)/15 }
  func[i+LEN_OTU,"Nleach"] <- new_nleach
  
  pleach_add <- round(rnorm(1,0,new_pleach_sd),digits=3)
  new_pleach <- func[i+LEN_OTU,"Pleach"] + pleach_add
  print(func[i+LEN_OTU,"Pleach"])
  print(pleach_add)
  print(new_pleach)
  if(new_pleach <= 0) { new_pleach <- rchisq(1,3)/15; print(new_pleach) }
  func[i+LEN_OTU,"Pleach"] <- new_pleach
  print("-------------")
  
  new_decom <- func[i+LEN_OTU,"Decom"] + rnorm(1,0,new_decom_sd)
  if(new_decom <= 0) { new_decom <- rchisq(1,3)/15 }
  func[i+LEN_OTU,"Decom"] <- new_decom
}


##### Change to model-matrix format #####

V <- MATLEN * (MATLEN + 1)/2

NEW.DATA <- matrix(0,dim(OTU_DATAT_RES)[1]+1,V)

CNM <- colnames(B.META)

mean.bact40 <- mean(Bact40)
sd.bact40 <- sd(Bact40)
sd.b.meta <- sd(B.META[as.matrix(B.META) != 0])

m <- 1
for(i in 2:(dim(OTU_DATAT_RES)[1]+1))
{
  l <- 1
  n <- 1
  for(j in 1:(dim(OTU_DATAT_RES)[2]))
  {
    for(k in j:dim(OTU_DATAT_RES)[2])
    {
      nmj <- CNM[j]
      nmk <- CNM[k]
      if (i == 2)
      {
        NEW.DATA[1,l] <- paste(nmj,"x",nmk)
      }
      
      if (nmj == nmk)
      {
        ### adjust relative abundances so they have mean 0 and same SD as other entries in B matrix
        ### don't use relative abundances, use presence/absence
        NEW.DATA[i,l] <-  ifelse(OTU_DATAT_RES[i-1,nmj] > 0, 1,
                                 #sd.b.meta*(Bact40[m,nmk] - mean.bact40)/sd.bact40,
                                 0)
      }
      else{
        NEW.DATA[i,l] <- ifelse((OTU_DATAT_RES[i-1,nmj] > 0) & (OTU_DATAT_RES[i-1,nmk] > 0), 
                               B.META[nmj,nmk],
                               0)
      }
      l <- l + 1
    }
  }
  m <- m+1
  if(i == 51)
  {
    m <- 1
  }
}

rownames(NEW.DATA) <- c("Edge",rownames(OTU_DATAT_RES))

XY <- cbind(NEW.DATA,c("Nleach",func$Nleach),c("Pleach",func$Pleach),c("Decom",func$Decom))

colnames(XY) <- XY[1,]
XY <- XY[-c(1),]
XY[1:10,180:200]

write.csv(XY,"Wagg_etal_data/cutoff=40_out=XY_100samples.csv",row.names = FALSE)


# Stats about networks

MATLEN <- dim(META_MATRIX)[1]
MATLEN

OTU_DATAT_RES.df <- data.frame(OTU_DATAT_RES)
OTU_DATAT_RES.df$num_nonzero <- rowSums(OTU_DATAT_RES.df!=0)
OTU_DATAT_RES.df$num_nonzero[1:50]
OTU_DATAT_RES.df$num_nonzero[51:100]

# ORIGINAL DATA
min(OTU_DATAT_RES.df$num_nonzero[1:50])
max(OTU_DATAT_RES.df$num_nonzero[1:50])
mean(OTU_DATAT_RES.df$num_nonzero[1:50])

# GENERATED DATA
min(OTU_DATAT_RES.df$num_nonzero[51:100])
max(OTU_DATAT_RES.df$num_nonzero[51:100])
mean(OTU_DATAT_RES.df$num_nonzero[51:100])

# Meta-Matrix
sum(upperTriangle(B.META,byrow=T) != 0) / length(upperTriangle(B.META,byrow=T))


### Stats about Pleach 
print("100 total samples")
XYs <- data.frame(XY)
XYs$Pleach <- as.numeric(XYs$Pleach)
ggplot(data=XYs) + geom_histogram(aes(x=Pleach))
ggplot(data=XYs[1:50,]) + geom_histogram(aes(x=Pleach))
ggplot(data=XYs[51:100,]) + geom_histogram(aes(x=Pleach))

summary(XYs$Pleach[1:50])
round(summary(XYs$Pleach[51:100]),4)


### 

mean(B.META[as.matrix(B.META) != 0])
min(B.META[as.matrix(B.META) != 0])
max(B.META[as.matrix(B.META) != 0])
sd(B.META[as.matrix(B.META) != 0])


## 40 cutoff - 200 total samples
set.seed(122522)
META_MATRIX <- read.csv("Wagg_etal_data/Meta_network_weighted_40.csv")

OTU_DATA <- read.csv("Wagg_etal_data/Bacterial_OTU_data.csv")
OTU_DATA[1:10,1:10]


OTU_DATAT <- transpose(OTU_DATA)
OTU_DATAT <- OTU_DATAT[-c(1),]

MATLEN <- dim(META_MATRIX)[1]

rownames(OTU_DATAT) <- colnames(OTU_DATA[,-c(1)])
colnames(OTU_DATAT) <- parse_number(OTU_DATA$X)

OTU_DATAT <- OTU_DATAT[,order(as.numeric(colnames(OTU_DATAT)))]
OTU_DATAT[1:10,1:10]

B.META <- META_MATRIX
dim(B.META)
B.META[1:10,1:10]

rownames(B.META) <- parse_number(gsub('\\.','',B.META$X))
B.META <- B.META[,c(-1)]
dim(B.META)
colnames(B.META) <- parse_number(gsub('\\.','',colnames(B.META)))
B.META[MATLEN-10:MATLEN,MATLEN-10:MATLEN]

B.META <- B.META[order(as.numeric(rownames(B.META))),order(as.numeric(colnames(B.META)))]
B.META[1:10,1:10]

upperTriangle(B.META,byrow=T)

colnames(Bact40) <- parse_number(gsub('\\.','',colnames(Bact40)))

OTU_DATAT_RES_TMP <- OTU_DATAT[,colnames(OTU_DATAT) %in% colnames(B.META)]
dim(OTU_DATAT_RES_TMP)

LEN_OTU <- dim(OTU_DATAT_RES_TMP)[1]
OTU_DATAT_RES <- data.frame(matrix(0,LEN_OTU*4,dim(OTU_DATAT_RES_TMP)[2]))
colnames(OTU_DATAT_RES) <- colnames(OTU_DATAT_RES_TMP)

func <- data.frame(matrix(0,LEN_OTU*2,dim(Func)[2]))
colnames(func) <- colnames(Func)

keep_prob <- 0.90

new_nleach_sd <- (sd(Func$Nleach[1:50]))^(1/50)
new_pleach_sd <- (sd(Func$Pleach[1:50]))/2
new_decom_sd <- (sd(Func$Decom[1:50]))^(1/50)

for(i in 1:LEN_OTU)
{
  num_necesary <- length(OTU_DATAT_RES_TMP[i,])
  OTU_DATAT_RES[i,] <- OTU_DATAT_RES_TMP[i,]
  for (j in 1:num_necesary)
  {
    OTU_DATAT_RES[i+LEN_OTU,j] <- ifelse(OTU_DATAT_RES_TMP[i,j] > 0, rbinom(1,1,keep_prob), rbinom(1,1,1-keep_prob))
    OTU_DATAT_RES[i+(2*LEN_OTU),j] <- ifelse(OTU_DATAT_RES_TMP[i,j] > 0, rbinom(1,1,keep_prob), rbinom(1,1,1-keep_prob))
    OTU_DATAT_RES[i+(3*LEN_OTU),j] <- ifelse(OTU_DATAT_RES_TMP[i,j] > 0, rbinom(1,1,keep_prob), rbinom(1,1,1-keep_prob))
  }
  func[i,] <- Func[i,]
  func[i+LEN_OTU,] <- Func[i,]
  func[i+(2*LEN_OTU),] <- Func[i,]
  func[i+(3*LEN_OTU),] <- Func[i,]
  
  new_nleach <- func[i+LEN_OTU,"Nleach"] + round(rnorm(1,0,new_nleach_sd),digits=2)
  if(new_nleach <= 0) { new_nleach <- rchisq(1,3)/15 }
  func[i+LEN_OTU,"Nleach"] <- new_nleach
  func[i+(2*LEN_OTU),"Nleach"] <- new_nleach
  func[i+(3*LEN_OTU),"Nleach"] <- new_nleach
  
  pleach_add <- round(rnorm(1,0,new_pleach_sd),digits=3)
  new_pleach <- func[i+LEN_OTU,"Pleach"] + pleach_add
  if(new_pleach <= 0) { new_pleach <- rchisq(1,3)/15; print(new_pleach) }
  func[i+LEN_OTU,"Pleach"] <- new_pleach
  func[i+(2*LEN_OTU),"Pleach"] <- new_pleach
  func[i+(3*LEN_OTU),"Pleach"] <- new_pleach
  
  
  new_decom <- func[i+LEN_OTU,"Decom"] + rnorm(1,0,new_decom_sd)
  if(new_decom <= 0) { new_decom <- rchisq(1,3)/15 }
  func[i+LEN_OTU,"Decom"] <- new_decom
  func[i+(2*LEN_OTU),"Decom"] <- new_decom
  func[i+(3*LEN_OTU),"Decom"] <- new_decom
}


##### Change to model-matrix format #####

V <- MATLEN * (MATLEN + 1)/2

NEW.DATA <- matrix(0,dim(OTU_DATAT_RES)[1]+1,V)

CNM <- colnames(B.META)

mean.bact40 <- mean(Bact40)
sd.bact40 <- sd(Bact40)
sd.b.meta <- sd(B.META[as.matrix(B.META) != 0])

m <- 1
for(i in 2:(dim(OTU_DATAT_RES)[1]+1))
{
  l <- 1
  n <- 1
  for(j in 1:(dim(OTU_DATAT_RES)[2]))
  {
    for(k in j:dim(OTU_DATAT_RES)[2])
    {
      nmj <- CNM[j]
      nmk <- CNM[k]
      if (i == 2)
      {
        NEW.DATA[1,l] <- paste(nmj,"x",nmk)
      }
      
      if (nmj == nmk)
      {
        ### adjust relative abundances so they have mean 0 and same SD as other entries in B matrix
        ### don't use relative abundances, use presence/absence
        NEW.DATA[i,l] <-  ifelse(OTU_DATAT_RES[i-1,nmj] > 0, 1,
                                 #sd.b.meta*(Bact40[m,nmk] - mean.bact40)/sd.bact40,
                                 0)
      }
      else{
        NEW.DATA[i,l] <- ifelse((OTU_DATAT_RES[i-1,nmj] > 0) & (OTU_DATAT_RES[i-1,nmk] > 0), 
                                B.META[nmj,nmk],
                                0)
      }
      l <- l + 1
    }
  }
  m <- m+1
  if(i == 51)
  {
    m <- 1
  }
}

rownames(NEW.DATA) <- c("Edge",rownames(OTU_DATAT_RES))

XY <- cbind(NEW.DATA,c("Nleach",func$Nleach),c("Pleach",func$Pleach),c("Decom",func$Decom))

colnames(XY) <- XY[1,]
XY <- XY[-c(1),]
XY[1:10,180:200]

write.csv(XY,"Wagg_etal_data/cutoff=40_out=XY_200samples.csv",row.names = FALSE)


# Stats about networks

MATLEN <- dim(META_MATRIX)[1]
MATLEN

OTU_DATAT_RES.df <- data.frame(OTU_DATAT_RES)
OTU_DATAT_RES.df$num_nonzero <- rowSums(OTU_DATAT_RES.df!=0)
OTU_DATAT_RES.df$num_nonzero[1:50]
OTU_DATAT_RES.df$num_nonzero[51:200]

# ORIGINAL DATA
min(OTU_DATAT_RES.df$num_nonzero[1:50])
max(OTU_DATAT_RES.df$num_nonzero[1:50])
mean(OTU_DATAT_RES.df$num_nonzero[1:50])

# GENERATED DATA
min(OTU_DATAT_RES.df$num_nonzero[51:200])
max(OTU_DATAT_RES.df$num_nonzero[51:200])
mean(OTU_DATAT_RES.df$num_nonzero[51:200])

# Meta-Matrix
sum(upperTriangle(B.META,byrow=T) != 0) / length(upperTriangle(B.META,byrow=T))


### Stats about Pleach 
print("200 total samples")
XYs <- data.frame(XY)
XYs$Pleach <- as.numeric(XYs$Pleach)
ggplot(data=XYs) + geom_histogram(aes(x=Pleach))
ggplot(data=XYs[1:50,]) + geom_histogram(aes(x=Pleach))
ggplot(data=XYs[51:200,]) + geom_histogram(aes(x=Pleach))

summary(XYs$Pleach[1:50])
round(summary(XYs$Pleach[51:200]),4)


