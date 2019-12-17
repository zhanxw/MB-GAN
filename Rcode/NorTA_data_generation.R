# note: generate simulated data by Normal-to-Anything (NorTA) #

# libaray to implement NorTA #
library(SpiecEasi)

# load the original count data #
load("NielsenHB_original_data.Rdata", verbose = T)

# input of NorTA should be a sample-by-taxon count matrix #
Y.species = t(NielsenHB_count)
dim(Y.species)

# separate for the case and control groups #
Y.count.case = Y.species[patient.info[["disease"]] == "IBD", ]
Y.count.ctrl = Y.species[patient.info[["disease"]] =="healthy", ]

#### (A) NorTA simulate samples for the case group ####
# need to first remove the sparse taxa that only have zeros across all the samples #
zero.count.by.taxa = apply(Y.count.case, 2, function(x){sum(x == 0)})
# number of taxa to be removed #
sum(zero.count.by.taxa==nrow(Y.count.case))

Y.count.case.filter = Y.count.case[, zero.count.by.taxa!=nrow(Y.count.case) ]

# simulate count data by NorTA #
nsimu = 1000
set.seed(2174171)
Y.NorTA.case = synth_comm_from_counts(comm = Y.count.case.filter, 
                                      mar=2, distr='zinegbin', n=nsimu)
colnames(Y.NorTA.case) = colnames(Y.count.case.filter)
dim(Y.NorTA.case)

# extract the species level count matrix #
case.count.NorTA = Y.NorTA.case
case.count.NorTA.spec = case.count.NorTA
case.count.NorTA.spec = case.count.NorTA.spec[, grep("s__", colnames(case.count.NorTA.spec))]

# convert to the compositional data #
case.compo.NorTA.spec = t(apply(case.count.NorTA.spec, 1, function(x){x/sum(x)}))


#### (B) NorTA simulate samples for the control group ####
# remove the too sparse taxa #
zero.count.by.taxa = apply(Y.count.ctrl, 2, function(x){sum(x == 0)})
sum(zero.count.by.taxa==nrow(Y.count.ctrl))

Y.count.ctrl.filter = Y.count.ctrl[, zero.count.by.taxa!=nrow(Y.count.ctrl) ]

# simulate count data by normal-to-anything #
Y.count.ctrl.tmp = Y.count.ctrl.filter
nsimu = 1000
set.seed(2174171)
Y.NorTA.ctrl = synth_comm_from_counts(comm = Y.count.ctrl.tmp, 
                                      mar=2, distr='zinegbin', n=nsimu)
colnames(Y.NorTA.ctrl) = colnames(Y.count.ctrl.filter)

# extract the species level count matrix #
ctrl.count.NorTA = Y.NorTA.ctrl
ctrl.count.NorTA.spec = ctrl.count.NorTA
ctrl.count.NorTA.spec = ctrl.count.NorTA.spec[, grep("s__", colnames(ctrl.count.NorTA.spec))]

# convert to the compositional data #
ctrl.compo.NorTA.spec = t(apply(ctrl.count.NorTA.spec, 1, function(x){x/sum(x)}))
dim(ctrl.compo.NorTA.spec)

#### comments ####
# The results are saved in the "data" folder 
# Can directly load the data "Data_by_NorTA.Rdata
