library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(gplots)
source("utility.R")
#
setwd("../data/")
# load the names by NorTA (this is a subset)
load("taxa_name_NorTA.Rdata")

# load the real data to find the top 60 most abundant taxa 
load("Data_NielsenHB_real.Rdata")

# 

abd.count = apply(case.real.spec, 2, function(x){sum(x!=0)})
abd.idx = names(rev(sort(abd.count))[1:60])

#### 1. case group ####
load("Data_by_NorTA.Rdata")

# check sparsity #
sparsity.case.NorTA = apply(case.compo.NorTA.spec, 1, get.sparsity)
summary(sparsity.case.NorTA); 
#round(var(sparsity.case.NorTA), 4)

# get the most abundant taxa (60)
# NorTA
NorTA.top.abd = case.compo.NorTA.spec[, abd.idx]
colnames(NorTA.top.abd) = NULL
rownames(NorTA.top.abd) = NULL
NorTA.top.df = melt(NorTA.top.abd)
colnames(NorTA.top.df) = c("row", "col", "Abundance")
ggplot(NorTA.top.df, aes(col,row, fill= Abundance)) + 
  geom_tile()+
  scale_fill_gradient(low="gray99", high="gray20",limits=c(0, 1))+
  xlab("Taxa") + ylab("Sample")+
  theme(panel.border = element_rect(colour = "grey20", fill=NA, size=1),
        panel.background = element_blank())+
  scale_x_continuous(limits = c(0, 60), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1000), expand = c(0, 0)) 


# real 
real.top.abd = case.real.spec[, abd.idx]
colnames(real.top.abd) = NULL
rownames(real.top.abd) = NULL
real.top.df = melt(real.top.abd)
colnames(real.top.df) = c("row", "col", "Abundance")
ggplot(real.top.df, aes(col,row, fill= Abundance)) + 
  geom_tile()+
  scale_fill_gradient(low="gray99", high="gray20",limits=c(0, 1))+
  xlab("Taxa") + ylab("Sample")+
  theme(panel.border = element_rect(colour = "grey20", fill=NA, size=1),
        panel.background = element_blank())+
  scale_x_continuous(limits = c(0, 60), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 148), expand = c(0, 0)) 

# MB-GAN
load("Data_by_WGAN_v1.Rdata")
MBGAN.top.abd = as.matrix(case.WGAN[, abd.idx])
colnames(MBGAN.top.abd) = NULL
rownames(MBGAN.top.abd) = NULL
MBGAN.top.df = melt(MBGAN.top.abd)
colnames(MBGAN.top.df) = c("row", "col", "Abundance")
ggplot(MBGAN.top.df, aes(col,row, fill= Abundance)) + 
  geom_tile()+
  scale_fill_gradient(low="gray99", high="gray20",limits=c(0, 1))+
  xlab("Taxa") + ylab("Sample")+
  theme(panel.border = element_rect(colour = "grey20", fill=NA, size=1),
        panel.background = element_blank())+
  scale_x_continuous(limits = c(0, 60), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1000), expand = c(0, 0)) 

