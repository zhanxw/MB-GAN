library(corrplot)
source("utility.R")

setwd("../data/")
# load the names by NorTA (this is a subset)
load("taxa_name_NorTA.Rdata")

#### 1. case group ####
# (a) load the real data #
load("Data_NielsenHB_real.Rdata")

# restrict to the subgroup by case group #
case.real.filter = case.real.spec[, colnames(case.real.spec) %in% taxa.case.NorTA]

# (b) WGAN #
load("Data_by_WGAN_v1.Rdata")

# restrict to the subgroup by case group #
case.WGAN.filter = case.WGAN[, colnames(case.WGAN) %in% taxa.case.NorTA]
# removed those extremely sparse taxa #
sparsity.case.WGAN.ori = apply(case.WGAN.filter, 1, get.sparsity)
sparsity.cutoff = 0.95
WGAN.sample.rm.idx = which(sparsity.case.WGAN.ori > sparsity.cutoff)
case.WGAN.filter = case.WGAN.filter[-WGAN.sample.rm.idx, ]

vis.spearman.mat(Xmat1 = case.real.filter, 
                 Xmat2 = case.WGAN.filter,
                 threshold = 0.1)

# (c) NorTA #
load("Data_by_NorTA.Rdata")
vis.spearman.mat(Xmat1 = case.real.filter,
                 Xmat2 = case.compo.NorTA.spec,
                 threshold = 0.1)

# distribution #
par(mar = c(4, 4, 1, 1))
par(mfrow = c(1, 2))
compare.spearman.hist (Xmat1 = case.real.filter,
                        Xmat2 = case.WGAN.filter,
                        Xmat3 = case.compo.NorTA.spec,
                        threshold = 0.1)



#### 2. contrl group ####
ctrl.real.filter = ctrl.real.spec[, colnames(ctrl.real.spec) %in% taxa.ctrl.NorTA]
# (a) WGAN vs TRUE #
ctrl.WGAN.filter = ctrl.WGAN[, colnames(ctrl.WGAN) %in% taxa.ctrl.NorTA]
sparsity.ctrl.WGAN.ori = apply(ctrl.WGAN.filter , 1, get.sparsity)
sparsity.cutoff = 0.95
WGAN.sample.rm.idx = which(sparsity.ctrl.WGAN.ori > sparsity.cutoff)
ctrl.WGAN.filter = ctrl.WGAN.filter[-WGAN.sample.rm.idx, ]


vis.spearman.mat(Xmat1 = ctrl.real.filter,
                 Xmat2 = ctrl.WGAN.filter,
                 threshold = 0.1)


# (b) NorTA #
vis.spearman.mat(Xmat1 = ctrl.real.filter,
                 Xmat2 = ctrl.compo.NorTA.spec,
                 threshold = 0.1)


# compare distributions #
par(mar = c(4, 4, 1, 1))
par(mfrow = c(1, 2))
compare.spearman.hist (Xmat1 = ctrl.real.filter,
                        Xmat2 = ctrl.WGAN.filter,
                        Xmat3 = ctrl.compo.NorTA.spec,
                        threshold = 0.1)


