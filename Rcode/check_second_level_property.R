source("utility.R")
library(corrplot)
# load the names by NorTA (this is a subset)
load("../data/taxa_name_NorTA.Rdata", verbose = T)

#### 1. case group ####
# (a) load the real data #
load("../data/Data_NielsenHB_real.Rdata", verbose = T)
load("../data/Data_by_WGAN_v1.Rdata", verbose = T)
load("../data/Data_by_NorTA.Rdata", verbose = T)
# new data #
load("../data/Data_by_metap_full.Rdata")

# make sure the abundance matrices have the same spcies #
case.real.filter = case.real.spec[, colnames(case.real.spec) %in% taxa.case.NorTA]
case.MB.filter = case.WGAN[, colnames(case.WGAN) %in% taxa.case.NorTA]
case.norta.filter = case.compo.NorTA.spec
case.metap.filter = case.metap.full[, colnames(case.metap.full)%in% taxa.case.NorTA]
set.seed(2147213)
case.metap.filter = case.metap.filter[sample(1:nrow(case.metap.filter),size = 1000), ]

# extract the taxa to be included in the plot #
# 1. remove the taxa with all zeros across all the samples #
zero.real = apply(case.real.filter, 2, function(x){sum(x == 0)}) == nrow(case.real.filter)
zero.MB = apply(case.MB.filter, 2, function(x){sum(x == 0)}) == nrow(case.MB.filter)
zero.norta = apply(case.norta.filter, 2, function(x){sum(x == 0)}) == nrow(case.norta.filter)
zero.metap = apply(case.metap.filter, 2, function(x){sum(x == 0)}) == nrow(case.metap.filter)

# 2. top abundant taxa in the real data #
threshold = 0.1
zero.count.by.taxa = apply(case.real.filter, 2, function(x){sum(x == 0)})
extr.idx = zero.count.by.taxa < (nrow(case.real.filter) * threshold)

# 3. get the name of taxa to be compared #
zero.rm.idx = zero.real | zero.MB |zero.norta | zero.metap
keep.idx = extr.idx & !zero.rm.idx
keep.nam = colnames(case.real.filter)[keep.idx]

# spearman correlation #
vis.spearman(Xmat1 = case.real.filter,
             Xmat2 = case.MB.filter,
             Xmat3 = case.norta.filter,
             Xmat4 = case.metap.filter,
             taxa.name = keep.nam, plot.hist = F,
             nbreak = 20)

scatter.spearman(Xmat1 = case.real.filter, 
                 Xmat2 = case.MB.filter, 
                 Xmat3 = case.norta.filter,
                 Xmat4 = case.metap.filter,
                 taxa.name = keep.nam)

# proportionality #
vis.proportionality(Xmat1 = case.real.filter,
                    Xmat2 = case.MB.filter,
                    Xmat3 = case.norta.filter,
                    Xmat4 = case.metap.filter,
                    taxa.name = keep.nam, nbreak = 20)

scatter.proportionality(Xmat1 = case.real.filter, 
                        Xmat2 = case.MB.filter, 
                        Xmat3 = case.norta.filter, 
                        Xmat4 = case.metap.filter,
                        taxa.name = keep.nam)


#### 2. control group ####

# make sure the abundance matrices have the same spcies #
ctrl.real.filter = ctrl.real.spec[, colnames(ctrl.real.spec) %in% taxa.ctrl.NorTA]
ctrl.MB.filter = ctrl.WGAN[, colnames(ctrl.WGAN) %in% taxa.ctrl.NorTA]
ctrl.norta.filter = ctrl.compo.NorTA.spec
ctrl.metap.filter = ctrl.metap.full[, colnames(ctrl.metap.full)%in% taxa.ctrl.NorTA]
set.seed(2147213)
ctrl.metap.filter = ctrl.metap.filter[sample(1:nrow(ctrl.metap.filter),size = 1000), ]

# extract the taxa to be included in the plot #
# 1. remove the taxa with all zeros across all the samples #
zero.real = apply(ctrl.real.filter, 2, function(x){sum(x == 0)}) == nrow(ctrl.real.filter)
zero.MB = apply(ctrl.MB.filter, 2, function(x){sum(x == 0)}) == nrow(ctrl.MB.filter)
zero.norta = apply(ctrl.norta.filter, 2, function(x){sum(x == 0)}) == nrow(ctrl.norta.filter)
zero.metap = apply(ctrl.metap.filter, 2, function(x){sum(x == 0)}) == nrow(ctrl.metap.filter)


# 2. top abundant taxa in the real data #
threshold = 0.1
zero.count.by.taxa = apply(ctrl.real.filter, 2, function(x){sum(x == 0)})
extr.idx = zero.count.by.taxa < (nrow(ctrl.real.filter) * threshold)

# 3. get the name of taxa to be compared #
zero.rm.idx = zero.real | zero.MB |zero.norta | zero.metap
keep.idx = extr.idx & !zero.rm.idx
keep.nam = colnames(ctrl.real.filter)[keep.idx]

# make the plot #
vis.spearman(Xmat1 = ctrl.real.filter,
             Xmat2 = ctrl.MB.filter,
             Xmat3 = ctrl.norta.filter,
             Xmat4 = ctrl.metap.filter,
             taxa.name = keep.nam, 
             plot.hist = F,
             nbreak = 20)

scatter.spearman(Xmat1 = ctrl.real.filter, 
                 Xmat2 = ctrl.MB.filter, 
                 Xmat3 = ctrl.norta.filter, 
                 Xmat4 = ctrl.metap.filter, 
                 taxa.name = keep.nam)

vis.proportionality(ctrl.real.filter,
                    ctrl.MB.filter,
                    ctrl.norta.filter,
                    ctrl.metap.filter, taxa.name = keep.nam, nbreak = 20)

scatter.proportionality(Xmat1 = ctrl.real.filter, 
                        Xmat2 = ctrl.MB.filter, 
                        Xmat3 = ctrl.norta.filter, 
                        Xmat4 = ctrl.metap.filter, 
                        taxa.name = keep.nam)
