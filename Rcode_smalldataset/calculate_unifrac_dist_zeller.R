library("phyloseq")
library(vegan)
library(GUniFrac)

get.unifrac <- function(obj.phylo) {
  tr <- phy_tree(obj.phylo)
  tmp.otu <- as.data.frame(t(otu_table(obj.phylo)@.Data))
  tmp.ret <- GUniFrac(tmp.otu, tr)
  cat("finish computing the unifrac distance \n")
  tmp.d = tmp.ret$unifracs[,,"d_UW"]
  metap.nmds = vegan::metaMDS(tmp.d, trymax = 500)
}

# load the names by NorTA (this is a subset)
load("../data/taxa_name_NorTA_zeller.Rdata", verbose = T)
load("../data/phylo_tree_info_zeller.Rdata", verbose = T)

#### 1. case group ####
# (a) load the real data #
load("../data/zeller_species_count.Rdata", verbose = T)
load("../data/MBGAN_zeller.Rdata", verbose = T)
load("../data/Data_by_NorTA_zeller.Rdata", verbose = T)
load("../data/Data_by_metap_zeller.Rdata", verbose = T)

# convert the real data into compositional formula #
case.real.spec = t(apply(Y_case, 1, function(x){x/sum(x)}))
ctrl.real.spec = t(apply(Y_ctrl, 1, function(x){x/sum(x)}))

# make sure the abundance matrices have the same spcies #
case.real.filter = case.real.spec[, colnames(case.real.spec) %in% taxa.case.NorTA]
case.MB.filter = abundance.case[, colnames(abundance.case) %in% taxa.case.NorTA]
case.norta.filter = case.compo.NorTA.spec
case.metap.filter = case.metap.full[, colnames(case.metap.full)%in% taxa.case.NorTA]
set.seed(2147213)
case.metap.filter = case.metap.filter[sample(1:nrow(case.metap.filter),size = 1000), ]

#### construct the phyloseq object ####

# real #
case.real.tab = otu_table(t(case.real.filter), taxa_are_rows = T)
case.real.phylo = phyloseq(case.real.tab, tree.tip)
nmds.case.real = get.unifrac(case.real.phylo)

# mb-gan #
case.MB.tab = otu_table(t(case.MB.filter), taxa_are_rows = T)
case.MB.phylo = phyloseq(case.MB.tab, tree.tip)
nmds.case.MB = get.unifrac(case.MB.phylo)

# metasparsim #
case.metap.tab = otu_table(t(case.metap.filter), taxa_are_rows = T)
case.metap.phylo = phyloseq(case.metap.tab, tree.tip)
nmds.case.metap = get.unifrac(case.metap.phylo)

# normal-to-anything # 
case.norta.tab = otu_table(t(case.norta.filter), taxa_are_rows = T)
case.norta.phylo = phyloseq(case.norta.tab, tree.tip)
nmds.case.norta = get.unifrac(case.norta.phylo)


#### 2. control group ####
# make sure the abundance matrices have the same spcies #
ctrl.real.filter = ctrl.real.spec[, colnames(ctrl.real.spec) %in% taxa.ctrl.NorTA]
ctrl.MB.filter = abundance.ctrl[, colnames(abundance.ctrl) %in% taxa.ctrl.NorTA]
ctrl.norta.filter = ctrl.compo.NorTA.spec
ctrl.metap.filter = ctrl.metap.full[, colnames(ctrl.metap.full)%in% taxa.ctrl.NorTA]
set.seed(2147213)
ctrl.metap.filter = ctrl.metap.filter[sample(1:nrow(ctrl.metap.filter),size = 1000), ]

#### construct the phyloseq object ####

# real #
ctrl.real.tab = otu_table(t(ctrl.real.filter), taxa_are_rows = T)
ctrl.real.phylo = phyloseq(ctrl.real.tab, tree.tip)
nmds.ctrl.real = get.unifrac(ctrl.real.phylo)

# mb-gan #
ctrl.MB.tab = otu_table(t(ctrl.MB.filter), taxa_are_rows = T)
ctrl.MB.phylo = phyloseq(ctrl.MB.tab, tree.tip)
nmds.ctrl.MB = get.unifrac(ctrl.MB.phylo)

# metasparsim #
ctrl.metap.tab = otu_table(t(ctrl.metap.filter), taxa_are_rows = T)
ctrl.metap.phylo = phyloseq(ctrl.metap.tab, tree.tip)
nmds.ctrl.metap = get.unifrac(ctrl.metap.phylo)

# normal-to-anything # 
ctrl.norta.tab = otu_table(t(ctrl.norta.filter), taxa_are_rows = T)
ctrl.norta.phylo = phyloseq(ctrl.norta.tab, tree.tip)
nmds.ctrl.norta = get.unifrac(ctrl.norta.phylo)

save.nam = paste0("../data/unifrac_nmds_zeller.Rdata")
save(nmds.ctrl.real, nmds.ctrl.MB, nmds.ctrl.metap, nmds.ctrl.norta, 
     nmds.case.real, nmds.case.MB, nmds.case.metap, nmds.case.norta, 
     file = save.nam)