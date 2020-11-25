source("utility.R")
load("../data/mgban_TestDifferentialAbundanceSimulation.Rdata")

#### differential abundance analysis on the labelded data ####
# match the name on the species level #
simu.specis.nam = colnames(combined.file)
simu.specis.nam = gsub("\\.", "|", simu.specis.nam)
real.taxa.nam.full = colnames(Y.case)
#real.taxa.nam.full = sub('.*\\|', '', real.taxa.nam.full)
species.match.idx = real.taxa.nam.full %in% simu.specis.nam
Y.case = Y.case[, species.match.idx]; 
Y.ctrl = Y.control[, species.match.idx]

## this is for checking if the simulated data can be used to conduct association analysis #
combined.mat = as.matrix(combined.file)
dim(combined.mat)
combined.mat[combined.mat<2e-4] = 0
pheno.idx.simu = rep(c(1, 0), each = 1000)
colnames(combined.mat) = simu.specis.nam

# create the matrix for the real data #
Y.real.combine = rbind(Y.case, Y.ctrl)
dim(Y.real.combine); range(Y.real.combine)
pheno.idx = rep(c(1, 0), c(nrow(Y.case), nrow(Y.ctrl)))

# remove the species with zero only #
zero.count.real = apply(Y.real.combine, 2, function(x){sum(x == 0)})
zero.rm.idx.real = zero.count.real == nrow(Y.real.combine)
zero.count.simu = apply(combined.mat, 2, function(x){sum(x == 0)})
zero.rm.idx.simu = zero.count.simu == nrow(combined.mat)

Y.simu.combine = combined.mat[, !(zero.rm.idx.simu | zero.rm.idx.real)]
Y.real.combine = Y.real.combine[, !(zero.rm.idx.simu | zero.rm.idx.real)]

# change the name to the bottom-most level #
colnames(Y.real.combine) = sub('.*\\|', '',colnames(Y.real.combine))
colnames(Y.simu.combine) = sub('.*\\|', '',colnames(Y.simu.combine))

# KW test #
KW.boxplot.top10 = KW.boxplot(Xmat1 = Y.real.combine, Xmat2 = Y.simu.combine, z1 = pheno.idx, z2 = pheno.idx.simu, 
                              top = 10, nresample = 100 )


