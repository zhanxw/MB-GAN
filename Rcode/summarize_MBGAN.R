setwd("../data/")
# summary data simulated by MB-GAN #
cutoff = 1e-4


# load the MB-GAN samples 
# 1) case group #
WGAN.case = read.table("stools_2_case_species.csv", 
                       sep = ',', header = T)
WGAN.case = WGAN.case[, -1]

# trim the small values #
WGAN.case[WGAN.case < cutoff] = 0


# 2) control group #
WGAN.ctrl = read.table("stools_2_ctrl_species.csv", 
                       sep = ',', header = T)
WGAN.ctrl = WGAN.ctrl[, -1]

# trim the small values #
WGAN.ctrl[WGAN.ctrl < cutoff] = 0


# rename all the taxa to the lowest rank #
name.case.ref = colnames(WGAN.case)
name.case.ref = sub('.*\\.', '',name.case.ref)
colnames(WGAN.case) = name.case.ref

name.ctrl.ref = colnames(WGAN.ctrl)
name.ctrl.ref = sub('.*\\.', '',name.ctrl.ref)
colnames(WGAN.ctrl) = name.ctrl.ref

rowSums(WGAN.case); rowSums(WGAN.ctrl)

# convert to compositional data #
case.WGAN = WGAN.case; ctrl.WGAN = WGAN.ctrl

# save the results #
save(case.WGAN , ctrl.WGAN, file = "Data_by_WGAN_v1.Rdata")
name.case.WGAN = name.case.ref; name.ctrl.WGAN = name.ctrl.ref
save(name.case.WGAN, name.ctrl.WGAN, file = "taxa_name_WGAN_v1.Rdata")

