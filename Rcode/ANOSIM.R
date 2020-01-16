library(phyloseq)
library(vegan)
setwd("../data/")

#### PERMANOVA case group ####
# filter the taxa based on NorTA results #
load("taxa_name_NorTA.Rdata")

# real data #
load("Data_NielsenHB_real.Rdata")
case.real.spec.filter = case.real.spec[, colnames(case.real.spec) %in% taxa.case.NorTA]

# MB-GAN #
load("Data_by_WGAN_v1.Rdata")
# restrict to the subgroup by case group #
case.WGAN.filter = case.WGAN[, colnames(case.WGAN) %in% taxa.case.NorTA]

# 1.1 permanova between real and WGAN (unfiltered) #
# check the order of taxa names #
all(colnames(case.real.spec) == colnames(case.WGAN))

case1.df = rbind(case.real.spec, case.WGAN)
case1.label = rep(c("real","WGAN"), c(nrow(case.real.spec), nrow(case.WGAN)))

# permanova #
permanova.case1 = adonis(case1.df ~ case1.label,permutations=1000, method = "bray")
print(as.data.frame(permanova.case1$aov.tab)["case1.label", "Pr(>F)"])
# [1] 0.1838162

# ANOSIM #
case.bray.1 = vegdist(x = case1.df, method="bray" )
case.ANOSIM.1 = anosim(case.bray.1, case1.label, permutations = 1000)
case.ANOSIM.1$signif
#0.3756244

# 1.2 permanova between real and WGAN (filtered) #
# check the order of taxa names #
all(colnames(case.real.spec.filter) == colnames(case.WGAN.filter))

case1f.df = rbind(case.real.spec.filter, case.WGAN.filter)
case1f.label = rep(c("real","WGAN"), c(nrow(case.real.spec.filter), nrow(case.WGAN.filter)))

permanova.case1f = adonis(case1f.df ~ case1f.label,permutations=1000, method = "bray")
print(as.data.frame(permanova.case1f$aov.tab)["case1f.label", "Pr(>F)"])
#[1] 0.2037962


# 1.3 permanova between real and NorTA #
load("Data_by_NorTA.Rdata")
all(colnames(case.real.spec.filter) == colnames(case.compo.NorTA.spec))

case2.df = rbind(case.real.spec.filter, case.compo.NorTA.spec)
case2.label = rep(c("real","NorTA"), c(nrow(case.real.spec.filter), nrow(case.compo.NorTA.spec)))


# permanova #
permanova.case2 = adonis(case2.df ~ case2.label,permutations=1000, method = "bray")
print(as.data.frame(permanova.case2$aov.tab)["case2.label", "Pr(>F)"])
#[1] 0.000999001

# ANOSIM #
case.bray.2 = vegdist(x = case2.df, method="bray" )
case.ANOSIM.2 = anosim(case.bray.2, case2.label, permutations = 1000)
case.ANOSIM.2$signif
#[1] 0.000999001



#### PERMANOVA control group ####

# 2.1 permanova between real and WGAN (unfiltered) #
all(colnames(ctrl.real.spec) == colnames(ctrl.WGAN))

ctrl1.df = rbind(ctrl.real.spec, ctrl.WGAN)
ctrl1.label = rep(c("real","WGAN"), c(nrow(ctrl.real.spec), nrow(ctrl.WGAN)))

# permanova #
permanova.ctrl1 = adonis(ctrl1.df ~ ctrl1.label,permutations=1000, method = "bray")
print(as.data.frame(permanova.ctrl1$aov.tab)["ctrl1.label", "Pr(>F)"])
# [1] 0.001998002

# ANOSIM #
ctrl.bray.1 = vegdist(x = ctrl1.df, method="bray" )
ctrl.ANOSIM.1 = anosim(ctrl.bray.1, ctrl1.label, permutations = 1000)
ctrl.ANOSIM.1$signif
#[1] 0.3006993

# 2.2 permanova between real and WGAN (filtered) #
ctrl.real.spec.filter = ctrl.real.spec[, colnames(ctrl.real.spec) %in% taxa.ctrl.NorTA]
ctrl.WGAN.filter = ctrl.WGAN[, colnames(ctrl.WGAN) %in% taxa.ctrl.NorTA]

all(colnames(ctrl.real.spec.filter) == colnames(ctrl.WGAN.filter))

ctrl1f.df = rbind(ctrl.real.spec.filter, ctrl.WGAN.filter)
ctrl1f.label = rep(c("real","WGAN"), c(nrow(ctrl.real.spec.filter), nrow(ctrl.WGAN.filter)))
permanova.ctrl1f = adonis(ctrl1f.df ~ ctrl1f.label,permutations=1000, method = "bray")
print(as.data.frame(permanova.ctrl1f$aov.tab)["ctrl1f.label", "Pr(>F)"])
#[1] 0.000999001


# 2.3 permanova between real and NorTA #
all(colnames(ctrl.real.spec.filter) == colnames(ctrl.compo.NorTA.spec))

ctrl2.df = rbind(ctrl.real.spec.filter, ctrl.compo.NorTA.spec)
ctrl2.label = rep(c("real","NorTA"), c(nrow(ctrl.real.spec.filter), nrow(ctrl.compo.NorTA.spec)))


# permanova #
permanova.ctrl2 = adonis(ctrl2.df ~ ctrl2.label,permutations=1000, method = "bray")
print(as.data.frame(permanova.ctrl2$aov.tab)["ctrl2.label", "Pr(>F)"])
#[1] 0.000999001

# ANOSIM #
ctrl.bray.2 = vegdist(x = ctrl2.df, method="bray" )
ctrl.ANOSIM.2 = anosim(ctrl.bray.2, ctrl2.label, permutations = 1000)
ctrl.ANOSIM.2$signif
#[1] 0.000999001



