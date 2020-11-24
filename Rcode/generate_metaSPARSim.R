# note: generate simulated data by metaSPARSimn  #
library(metaSPARSim)
library(GMPR)

load("../data/NielsenHB_otu.Rdata", verbose = T)
load("../data/NielsenHB_OTUinfo.Rdata", verbose = T)
load("../data/Data_NielsenHB_real.Rdata", verbose = T )


otu_species = otu_full[rownames(otu_full) %in% colnames(case.real.spec), ]

condition.list = vector("list", 2)
condition.list[[1]] = which(group.label == 1)
condition.list[[2]] = which(group.label == 2)
names(condition.list) = c("control", "IBD")

otu_species1 = otu_species[, group.label == 1]
condition.list1 = list(condition.list[["control"]])
names(condition.list1) = "control"

# gmpr normalized data #
gmpr.size.factor = GMPR(t(otu_species1/10000))
gmpr.norm.data = otu_species1
for(jj in 1:ncol(otu_species1)){
  gmpr.norm.data[, jj] = gmpr.norm.data[, jj] / gmpr.size.factor[jj]
}

metaSPARSim.para = estimate_parameter_from_data( round(otu_species1/10000), gmpr.norm.data, condition.list1)
names(metaSPARSim.para) = names(condition.list)

# simulate data #
sim_data = metaSPARSim(metaSPARSim.para)

# process the count matrix and get the compositional data #
data.metap = sim_data$counts
data.real = rbind(ctrl.real.spec, case.real.spec)

# 1. case group #
case.metap = t(data.metap[, grep("IBD", colnames(data.metap))])

# convert to compositional data #
case.metap = t(apply(case.metap, 1, function(x){x/sum(x)}))
reorder.idx = order(match(x = colnames(case.metap), table = colnames(case.real.spec)))
case.metap = case.metap[, reorder.idx]

# 2. ctrl group #
ctrl.metap = t(data.metap[, grep("control", colnames(data.metap))])

# convert to compositional data #
ctrl.metap = t(apply(ctrl.metap, 1, function(x){x/sum(x)}))
reorder.idx = order(match(x = colnames(ctrl.metap), table = colnames(ctrl.real.spec)))
ctrl.metap = ctrl.metap[, reorder.idx]

