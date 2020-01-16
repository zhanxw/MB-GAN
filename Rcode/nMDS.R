library(phyloseq)
library(vegan)
setwd("../data/")

# filter the taxa based on NorTA results #
load("taxa_name_NorTA.Rdata")
# real data #
load("Data_NielsenHB_real.Rdata")
case.real.spec.filter = case.real.spec[, colnames(case.real.spec) %in% taxa.case.NorTA]

#### case group with bray curtis distance####

#### 1.1 real data ####
real.nmds = metaMDS(case.real.spec.filter, distance="bray", k=2, trymax=1000)

#### 1.2 MB-GAN ####
load("Data_by_WGAN_v1.Rdata")
# restrict to the subgroup by case group #
case.WGAN.filter = case.WGAN[, colnames(case.WGAN) %in% taxa.case.NorTA]

wgan.nmds = metaMDS(case.WGAN.filter, distance="bray", k=2, trymax=500)
#points(wgan.nmds, display="sites", pch=20, col = alpha("red", 0.4))

#### 1.3 normal-to-anything ####
load("Data_by_NorTA.Rdata")
NorTA.nmds = metaMDS(case.compo.NorTA.spec, distance="bray", k=2, trymax = 500)

#### control group with bray curtis distance####
ctrl.real.spec.filter = ctrl.real.spec[, colnames(ctrl.real.spec) %in% taxa.ctrl.NorTA]

#### 2.1 real data ####
real.nmds.ctrl = metaMDS(ctrl.real.spec.filter, distance="bray", k=2, trymax=1000)

#### 2.2 MB-GAN ####
# restrict to the subgroup by ctrl group #
ctrl.WGAN.filter = ctrl.WGAN[, colnames(ctrl.WGAN) %in% taxa.ctrl.NorTA]

wgan.nmds.ctrl = metaMDS(ctrl.WGAN.filter, distance="bray", k=2, trymax = 500)

#### 2.3 normal-to-anything ####
NorTA.nmds.ctrl = metaMDS(ctrl.compo.NorTA.spec, distance="bray", k=2, trymax=500)

# save the output #
save(real.nmds, wgan.nmds,NorTA.nmds,
     real.nmds.ctrl,wgan.nmds.ctrl,NorTA.nmds.ctrl, 
     file = "MDSoutput_update_v1")























