source("../Rcode/utility.R")
library(cowplot)
library(ggpubr)
theme_set(theme_cowplot())
library(phyloseq)
library(vegan)

load("../data/taxa_name_NorTA_zeller.Rdata")
my_comparisons = list( c("Real", "MB-GAN"), 
                       c("Real", "NorTA"),  
                       c("Real", "metaSPARSim"))
#### 1. case group ####

#### 1.1 (a) sample sparsity ####
# (a) load the real data #
load("zeller_species_count.Rdata", verbose = T)
load("MBGAN_zeller.Rdata", verbose = T)
load("Data_by_NorTA_zeller.Rdata", verbose = T)
load("Data_by_metap_zeller.Rdata", verbose = T)

# convert the real data into compositional formula #
case.real.spec = t(apply(Y_case, 1, function(x){x/sum(x)}))
ctrl.real.spec = t(apply(Y_ctrl, 1, function(x){x/sum(x)}))

# reorder the columns (taxa) of the mbgan samples #
abundance.case = abundance.case[, colnames(abundance.case) %in% colnames(case.real.spec)]
reorder.idx = match(x = colnames(abundance.case), table = colnames(case.real.spec))
abundance.case = abundance.case[, order(reorder.idx)]

abundance.ctrl = abundance.ctrl[, colnames(abundance.ctrl) %in% colnames(ctrl.real.spec)]
reorder.idx = match(x = colnames(abundance.ctrl), table = colnames(ctrl.real.spec))
abundance.ctrl = abundance.ctrl[, order(reorder.idx)]

# restrict to the subgroup by case group #
case.real.filter = case.real.spec[, colnames(case.real.spec) %in% taxa.case.NorTA]
sparsity.case.real = apply(case.real.filter, 1, get.sparsity)
summary(sparsity.case.real); 

# (b) MB-GAN #

# restrict to the subgroup by case group #
case.WGAN.filter = abundance.case[, colnames(abundance.case) %in% taxa.case.NorTA]
# removed those extremely sparse taxa #
sparsity.case.WGAN.ori = apply(case.WGAN.filter, 1, get.sparsity)
sparsity.cutoff = 0.95
WGAN.sample.rm.idx = which(sparsity.case.WGAN.ori > sparsity.cutoff)
case.WGAN.filter = case.WGAN.filter[-WGAN.sample.rm.idx, ]

# check sparsity #
sparsity.case.WGAN = apply(case.WGAN.filter, 1, get.sparsity)
summary(sparsity.case.WGAN); 

# (c) NorTA #

# check sparsity #
sparsity.case.NorTA = apply(case.compo.NorTA.spec, 1, get.sparsity)
summary(sparsity.case.NorTA); 

# (d) metaSPARSim #

# check sparsity #
case.metap.filter = case.metap.full[, colnames(case.metap.full)%in% taxa.case.NorTA]
set.seed(2147213)
case.metap.filter = case.metap.filter[sample(1:nrow(case.metap.filter),size = 1000), ]
sparsity.case.metap = apply(case.metap.filter, 1, get.sparsity)


##### compare the general abundances #####

abd.data.frame = rbind(case.real.filter, case.WGAN.filter, case.compo.NorTA.spec, case.metap.filter)
abd.data.frame = data.frame(abd.data.frame)

#### 1.1 (b) species sparsity ####
case.spar.real = apply(case.real.filter, 2, get.sparsity)

par(mfrow = c(1, 2) )
abd.box1  = compare.abundance.boxplot(sparsity.real = case.spar.real, 
                                      df1 = case.real.filter, 
                                      df2 = case.WGAN.filter, 
                                      df3 = case.compo.NorTA.spec, 
                                      df4 = case.metap.filter, 
                                      0, 0.1)

compare.abundance.boxplot(sparsity.real = case.spar.real, 
                          df1 = case.real.filter, 
                          df2 = case.WGAN.filter, 
                          df3 = case.compo.NorTA.spec, 
                          df4 = case.metap.filter, 
                          0.1, 0.2, y.hi = NULL)

abd.box2 = compare.abundance.boxplot(sparsity.real = case.spar.real, 
                                     df1 = case.real.filter, 
                                     df2 = case.WGAN.filter, 
                                     df3 = case.compo.NorTA.spec, 
                                     df4 = case.metap.filter, 
                                     0.1, 0.2)
plot_grid(abd.box1, abd.box2, 
          labels = c("(a)", "(b)"),
          ncol = 2, nrow = 1)

#### 1.2 shannon index ####
# (a) real data #
shannon.case.real = apply(case.real.filter, 1, get.shannon.entropy)

# (b) WGAN #
shannon.case.WGAN = apply(case.WGAN.filter, 1, get.shannon.entropy)

# (c) NorTA #
shannon.case.NorTA = apply(case.compo.NorTA.spec, 1, get.shannon.entropy)

# (d) metaSPARSim #
shannon.case.metap = apply(case.metap.filter, 1, get.shannon.entropy)

case.diversity.df = data.frame(Type = rep(c("Real", "MB-GAN", "NorTA", "metaSPARSim"), 
                                          c(length(shannon.case.real), length(shannon.case.WGAN), 1000, length(shannon.case.metap))),
                               Diversity = c(shannon.case.real, 
                                             shannon.case.WGAN, 
                                             shannon.case.NorTA, 
                                             shannon.case.metap), stringsAsFactors = F)
case.diversity.df$Type = factor(case.diversity.df$Type, levels = c("Real", "MB-GAN", "NorTA", "metaSPARSim")) 

case.diversity = ggplot(case.diversity.df, aes(x=Type, y=Diversity, color = Type)) + 
  geom_boxplot(outlier.size = 0.3) + 
  xlab(" ") + ylim(0, 7)+ylab("Shannon Index")+
  scale_color_manual(values=c("black", rgb(1,0,0,0.5), rgb(0,0.5,0,0.5), rgb(0.12, 0.56, 1,0.5)))+
  theme(legend.position="none",
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  stat_compare_means(comparisons = list( c("Real", "MB-GAN"), 
                                         c("Real", "NorTA"),  
                                         c("Real", "metaSPARSim")), label = "p.signif")
plot(case.diversity)

#### 2. contrl group ####

#### 2.1 (A) sample sparsity ####

# (a) real
ctrl.real.filter = ctrl.real.spec[, colnames(ctrl.real.spec) %in% taxa.ctrl.NorTA]
sparsity.ctrl.real = apply(ctrl.real.filter, 1, get.sparsity)
summary(sparsity.ctrl.real); 

# (b) WGAN #
ctrl.WGAN.filter = abundance.ctrl[, colnames(abundance.ctrl) %in% taxa.ctrl.NorTA]
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

# (d) metaSPARSim #
# check sparsity #
ctrl.metap.filter = ctrl.metap.full[, colnames(ctrl.metap.full)%in% taxa.ctrl.NorTA]
set.seed(2147213)
ctrl.metap.filter = ctrl.metap.filter[sample(1:nrow(ctrl.metap.filter),size = 1000), ]
sparsity.ctrl.metap = apply(ctrl.metap.filter, 1, get.sparsity)
summary(sparsity.ctrl.metap); 

#### 2.1 (B) species sparsity ####
ctrl.spar.real = apply(ctrl.real.filter, 2, get.sparsity)
abd.box1 = compare.abundance.boxplot(sparsity.real = ctrl.spar.real, 
                                     df1 = ctrl.real.filter, 
                                     df2 = ctrl.WGAN.filter, 
                                     df3 = ctrl.compo.NorTA.spec, 
                                     df4 = ctrl.metap.filter, 
                                     0, 0.1, 
                                     y.hi = NULL)

abd.box2 = compare.abundance.boxplot(sparsity.real = ctrl.spar.real, 
                                     df1 = ctrl.real.filter, 
                                     df2 = ctrl.WGAN.filter, 
                                     df3 = ctrl.compo.NorTA.spec, 
                                     df4 = ctrl.metap.filter, 
                                     0.1, 0.2, 
                                     y.hi = NULL)

plot_grid(abd.box1, abd.box2, 
          labels = c("(a)", "(b)"),
          ncol = 2, nrow = 1)


#### 2.2 shannon index ####
# (a) real data #
shannon.ctrl.real = apply(ctrl.real.filter, 1, get.shannon.entropy)

# (b) WGAN #
shannon.ctrl.WGAN = apply(ctrl.WGAN.filter, 1, get.shannon.entropy)

# (c) NorTA #
shannon.ctrl.NorTA = apply(ctrl.compo.NorTA.spec, 1, get.shannon.entropy)

# (d) metaSPARSim #
shannon.ctrl.metap = apply(ctrl.metap.filter, 1, get.shannon.entropy)

ctrl.diversity.df = data.frame(Type = rep(c("Real", "MB-GAN", "NorTA", "metaSPARSim"), 
                                          c(length(shannon.ctrl.real), length(shannon.ctrl.WGAN), 1000, length(shannon.ctrl.metap))),
                               Diversity = c(shannon.ctrl.real, 
                                             shannon.ctrl.WGAN, 
                                             shannon.ctrl.NorTA, 
                                             shannon.ctrl.metap), stringsAsFactors = F)
ctrl.diversity.df$Type = factor(ctrl.diversity.df$Type, levels = c("Real", "MB-GAN", "NorTA", "metaSPARSim")) 

ctrl.diversity = ggplot(ctrl.diversity.df, aes(x=Type, y=Diversity, color = Type)) + 
  geom_boxplot(outlier.size = 0.3) + 
  xlab(" ") + ylim(0, 7)+ylab("Shannon Index")+
  scale_color_manual(values=c("black", rgb(1,0,0,0.5), rgb(0,0.5,0,0.5), rgb(0.12, 0.56, 1,0.5)))+
  theme(legend.position="none",
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  stat_compare_means(comparisons = list( c("Real", "MB-GAN"), 
                                         c("Real", "NorTA"),  
                                         c("Real", "metaSPARSim")), label = "p.signif")
plot(ctrl.diversity)


##### 3. nMDS #####
load("../data/unifrac_nmds_zeller.Rdata", verbose = T)

#### case group with unifrac distance####
ord.real = nmds.case.real[["points"]]; 
ord.MB = nmds.case.MB[["points"]]; 
ord.norta = nmds.case.norta[["points"]]; 
ord.metap = nmds.ctrl.metap[["points"]]; 


## case 1 ##
xrange = c(min(ord.real[, 1],ord.metap[, 1], ord.MB[, 1], ord.norta[, 1] ) - 0.01, 
           max(ord.real[, 1],ord.metap[, 1], ord.MB[, 1], ord.norta[, 1]  ) + 0.01)
yrange = c(min(ord.real[, 2],ord.metap[, 2], ord.MB[, 2], ord.norta[, 2]  ) - 0.01, 
           max(ord.real[, 2],ord.metap[, 2], ord.MB[, 2], ord.norta[, 2]  ) + 0.01)

par(mar = c(4, 6, 2, 1))
par(mfcol = c(1, 3))
plot(1, 1, type="n", main="Real vs MB-GAN",
     xlim = xrange, ylim = yrange, 
     xlab = "nMDS1", ylab = "nMDS2", cex.lab=2, cex.main=2)
points(ord.MB[, 1], ord.MB[, 2], pch=20, col = alpha(rgb(1,0,0,0.5), 0.5), cex = 0.8)
points(ord.real[, 1], ord.real[, 2], pch=20, col = "black")


plot(1, 1, type="n", main="Real vs Noraml-To-Anything (NorTA)",
     xlim = xrange, ylim = yrange, 
     xlab = "nMDS1", ylab = "nMDS2", cex.lab=2, cex.main=2)
points(ord.norta[, 1], ord.norta[, 2], pch=20, col = alpha(rgb(0,0.5,0,0.5), 0.5), cex = 0.8)
points(ord.real[, 1], ord.real[, 2], pch=20, col = "black")


plot(1, 1, type="n", main="Real vs MetaSPARSim",
     xlim = xrange, ylim = yrange, 
     xlab = "nMDS1", ylab = "nMDS2", cex.lab=2, cex.main=2)
points(ord.metap[, 1], ord.metap[, 2], pch=20, col = alpha(rgb(0.12, 0.56, 1,0.5), 0.5), cex = 0.8)
points(ord.real[, 1], ord.real[, 2], pch=20, col = "black")
# nmds_zeller_case: save size 15 * 5 #

#### control group with bray curtis distance####
ord.real = nmds.ctrl.real[["points"]]; 
ord.MB = nmds.ctrl.MB[["points"]]; 
ord.norta = nmds.ctrl.norta[["points"]]; 
ord.metap = nmds.ctrl.metap[["points"]]; 

xrange = c(min(ord.real[, 1],ord.metap[, 1], ord.MB[, 1], ord.norta[, 1] ) - 0.01, 
           max(ord.real[, 1],ord.metap[, 1], ord.MB[, 1], ord.norta[, 1]  ) + 0.01)
yrange = c(min(ord.real[, 2],ord.metap[, 2], ord.MB[, 2], ord.norta[, 2]  ) - 0.01, 
           max(ord.real[, 2],ord.metap[, 2], ord.MB[, 2], ord.norta[, 2]  ) + 0.01)

par(mar = c(4, 6, 2, 1))
par(mfcol = c(1, 3))
plot(1, 1, type="n", main="Real vs MB-GAN",
     xlim = xrange, ylim = yrange, 
     xlab = "nMDS1", ylab = "nMDS2", cex.lab=2, cex.main=2)
points(ord.MB[, 1], ord.MB[, 2], pch=20, col = alpha(rgb(1,0,0,0.5), 0.5), cex = 0.8)
points(ord.real[, 1], ord.real[, 2], pch=20, col = "black")


plot(1, 1, type="n", main="Real vs Noraml-To-Anything (NorTA)",
     xlim = xrange, ylim = yrange, 
     xlab = "nMDS1", ylab = "nMDS2", cex.lab=2, cex.main=2)
points(ord.norta[, 1], ord.norta[, 2], pch=20, col = alpha(rgb(0,0.5,0,0.5), 0.5), cex = 0.8)
points(ord.real[, 1], ord.real[, 2], pch=20, col = "black")

plot(1, 1, type="n", main="Real vs MetaSPARSim",
     xlim = xrange, ylim = yrange, 
     xlab = "nMDS1", ylab = "nMDS2", cex.lab=2, cex.main=2)
points(ord.metap[, 1], ord.metap[, 2], pch=20, col = alpha(rgb(0.12, 0.56, 1,0.5), 0.5), cex = 0.8)
points(ord.real[, 1], ord.real[, 2], pch=20, col = "black")

