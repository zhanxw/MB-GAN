library(cowplot)
library(ggpubr)
theme_set(theme_cowplot())
source("utility.R")

setwd("../data/")
# load the names by NorTA (this is a subset)
load("taxa_name_NorTA.Rdata")
my_comparisons = list( c("Real", "MB-GAN"), c("Real", "NorTA") )
#### 1. case group ####

#### 1.1 (a) sample sparsity ####
# (a) load the real data #
load("Data_NielsenHB_real.Rdata")

# restrict to the subgroup by case group #
case.real.filter = case.real.spec[, colnames(case.real.spec) %in% taxa.case.NorTA]
sparsity.case.real = apply(case.real.filter, 1, get.sparsity)
summary(sparsity.case.real); 

# (b) MB-GAN #
load("Data_by_WGAN_v1.Rdata")

# restrict to the subgroup by case group #
case.WGAN.filter = case.WGAN[, colnames(case.WGAN) %in% taxa.case.NorTA]
# removed those extremely sparse taxa #
sparsity.case.WGAN.ori = apply(case.WGAN.filter, 1, get.sparsity)
sparsity.cutoff = 0.95
WGAN.sample.rm.idx = which(sparsity.case.WGAN.ori > sparsity.cutoff)
case.WGAN.filter = case.WGAN.filter[-WGAN.sample.rm.idx, ]

# check sparsity #
sparsity.case.WGAN = apply(case.WGAN.filter, 1, get.sparsity)
summary(sparsity.case.WGAN); 

# (c) NorTA #
load("Data_by_NorTA.Rdata")
# check sparsity #
sparsity.case.NorTA = apply(case.compo.NorTA.spec, 1, get.sparsity)
summary(sparsity.case.NorTA); 

#### 1.1 (b) species sparsity ####
case.spar.real = apply(case.real.filter, 2, get.sparsity)

compare.abundance.boxplot(sparsity.real = case.spar.real, 
                          df1 = case.real.filter, 
                          df2 = case.WGAN.filter, 
                          df3 = case.compo.NorTA.spec, 
                          0, 0.1)


compare.abundance.boxplot(sparsity.real = case.spar.real, 
                          df1 = case.real.filter, 
                          df2 = case.WGAN.filter, 
                          df3 = case.compo.NorTA.spec, 
                          0.1, 0.2)

compare.abundance.boxplot(sparsity.real = case.spar.real, 
                          df1 = case.real.filter, 
                          df2 = case.WGAN.filter, 
                          df3 = case.compo.NorTA.spec, 
                          0.2, 0.3)

compare.abundance.boxplot(sparsity.real = case.spar.real, 
                          df1 = case.real.filter, 
                          df2 = case.WGAN.filter, 
                          df3 = case.compo.NorTA.spec, 
                          0.3, 0.4)

compare.abundance.boxplot(sparsity.real = case.spar.real, 
                          df1 = case.real.filter, 
                          df2 = case.WGAN.filter, 
                          df3 = case.compo.NorTA.spec, 
                          0.4, 0.5)

#### 1.2 shannon index ####
# (a) real data #
shannon.case.real = apply(case.real.filter, 1, get.shannon.entropy)

# (b) WGAN #
shannon.case.WGAN = apply(case.WGAN.filter, 1, get.shannon.entropy)

# (c) NorTA #
shannon.case.NorTA = apply(case.compo.NorTA.spec, 1, get.shannon.entropy)

case.diversity.df = data.frame(Type = rep(c("Real", "MBGAN", "NorTA"), c(length(shannon.case.real), length(shannon.case.WGAN), 1000)),
                               Diversity = c(shannon.case.real, shannon.case.WGAN, shannon.case.NorTA), stringsAsFactors = F)
case.diversity.df$Type = factor(case.diversity.df$Type, levels = c("Real", "MBGAN", "NorTA")) 

case.diversity = ggplot(case.diversity.df, aes(x=Type, y=Diversity, color = Type)) + 
  geom_boxplot(outlier.size = 0.5) + 
  xlab(" ") + ylim(0, 6)+ylab("Shannon Index")+
  scale_color_manual(values=c("blue", rgb(1,0,0,0.5), rgb(0,1,0,0.5)))+
  theme(legend.position="none",
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")
plot(case.diversity)


#### 2. contrl group ####

#### 2.1 (A) sample sparsity ####

# (a) real
ctrl.real.filter = ctrl.real.spec[, colnames(ctrl.real.spec) %in% taxa.ctrl.NorTA]
sparsity.ctrl.real = apply(ctrl.real.filter, 1, get.sparsity)
summary(sparsity.ctrl.real); 

# (b) WGAN #
ctrl.WGAN.filter = ctrl.WGAN[, colnames(ctrl.WGAN) %in% taxa.ctrl.NorTA]
# removed those extremely sparse taxa #
sparsity.ctrl.WGAN.ori = apply(ctrl.WGAN.filter , 1, get.sparsity)
sparsity.cutoff = 0.95
WGAN.sample.rm.idx = which(sparsity.ctrl.WGAN.ori > sparsity.cutoff)
ctrl.WGAN.filter = ctrl.WGAN.filter[-WGAN.sample.rm.idx, ]

sparsity.ctrl.WGAN = apply(ctrl.WGAN.filter, 1, get.sparsity)
summary(sparsity.ctrl.WGAN);


# (c) NorTA #
sparsity.ctrl.NorTA = apply(ctrl.compo.NorTA.spec, 1, get.sparsity)
summary(sparsity.ctrl.NorTA);

#### 2.1 (B) species sparsity ####
ctrl.spar.real = apply(ctrl.real.filter, 2, get.sparsity)
compare.abundance.boxplot(sparsity.real = ctrl.spar.real, 
                          df1 = ctrl.real.filter, 
                          df2 = ctrl.WGAN.filter, 
                          df3 = ctrl.compo.NorTA.spec, 
                          0, 0.1, 
                          y.hi = NULL)

compare.abundance.boxplot(sparsity.real = ctrl.spar.real, 
                          df1 = ctrl.real.filter, 
                          df2 = ctrl.WGAN.filter, 
                          df3 = ctrl.compo.NorTA.spec, 
                          0.1, 0.2)

compare.abundance.boxplot(sparsity.real = ctrl.spar.real, 
                          df1 = ctrl.real.filter, 
                          df2 = ctrl.WGAN.filter, 
                          df3 = ctrl.compo.NorTA.spec, 
                          0.2, 0.3)

compare.abundance.boxplot(sparsity.real = ctrl.spar.real, 
                          df1 = ctrl.real.filter, 
                          df2 = ctrl.WGAN.filter, 
                          df3 = ctrl.compo.NorTA.spec, 
                          0.3, 0.4)

compare.abundance.boxplot(sparsity.real = ctrl.spar.real, 
                          df1 = ctrl.real.filter, 
                          df2 = ctrl.WGAN.filter, 
                          df3 = ctrl.compo.NorTA.spec, 
                          0.4, 0.5)


#### 2.2 shannon index ####
# (a) real data #
shannon.ctrl.real = apply(ctrl.real.filter, 1, get.shannon.entropy)

# (b) WGAN #
shannon.ctrl.WGAN = apply(ctrl.WGAN.filter, 1, get.shannon.entropy)

# (c) NorTA #
shannon.ctrl.NorTA = apply(ctrl.compo.NorTA.spec, 1, get.shannon.entropy)

ctrl.diversity.df = data.frame(Type = rep(c("Real", "MBGAN", "NorTA"),
                                          c(length(shannon.ctrl.real), length(shannon.ctrl.WGAN), 1000)),
                               Diversity = c(shannon.ctrl.real, shannon.ctrl.WGAN, shannon.ctrl.NorTA), stringsAsFactors = F)
ctrl.diversity.df$Type = factor(ctrl.diversity.df$Type, levels = c("Real", "MBGAN", "NorTA")) 


ctrl.diversity = ggplot(ctrl.diversity.df, aes(x=Type, y=Diversity, color = Type)) + 
  geom_boxplot(outlier.size = 0.5) + xlab(" ") + ylim(0, 6)+ylab("Shannon Index")+
  scale_color_manual(values=c("blue", rgb(1,0,0,0.5), rgb(0,1,0,0.5)))+
  theme(legend.position="none",
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  stat_compare_means(comparisons = my_comparisons,  label = "p.signif")
plot(ctrl.diversity)




