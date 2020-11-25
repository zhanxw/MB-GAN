# note: generate simulated data by NorTA based on Zeller's study#

library(SpiecEasi)
source("../Rcode/utility.R")

#### compare with the performance of normal-to-anything (only performed on the case group ) ####
load("../data/zeller_species_count.Rdata")

Y.count.case = Y_case
Y.count.ctrl = Y_ctrl

#### (A) NorTA case group ####
# remove the too sparse taxa #
zero.count.by.taxa = apply(Y.count.case, 2, function(x){sum(x == 0)})
Y.count.case.filter = Y.count.case[, zero.count.by.taxa!=nrow(Y.count.case) ]

# simulate count data by normal-to-anything #
Y.count.case.tmp = round(Y.count.case.filter/10)
nsimu = 1000
set.seed(2174171)
Y.NorTA.case <- synth_comm_from_counts(comm = Y.count.case.tmp, 
                                       mar=2, distr='zinegbin', n=nsimu)
colnames(Y.NorTA.case) = colnames(Y.count.case.filter)

case.count.NorTA = Y.NorTA.case
case.count.NorTA.spec = case.count.NorTA
case.compo.NorTA.spec = t(apply(case.count.NorTA.spec, 1, function(x){x/sum(x)}))

#### (B) NorTA control group ####
# remove the too sparse taxa #
zero.count.by.taxa = apply(Y.count.ctrl, 2, function(x){sum(x == 0)})
Y.count.ctrl.filter = Y.count.ctrl[, zero.count.by.taxa!=nrow(Y.count.ctrl) ]

# simulate count data by normal-to-anything #
Y.count.ctrl.tmp = round(Y.count.ctrl.filter/10)
nsimu = 1000
set.seed(2174171)
Y.NorTA.ctrl <- synth_comm_from_counts(comm = Y.count.ctrl.tmp, 
                                       mar=2, distr='zinegbin', n=nsimu)
colnames(Y.NorTA.ctrl) = colnames(Y.count.ctrl.filter)

ctrl.count.NorTA = Y.NorTA.ctrl
ctrl.count.NorTA.spec = ctrl.count.NorTA
ctrl.compo.NorTA.spec = t(apply(ctrl.count.NorTA.spec, 1, function(x){x/sum(x)}))

save(Y.NorTA.case, case.count.NorTA, case.count.NorTA.spec, case.compo.NorTA.spec, 
     Y.NorTA.ctrl, ctrl.count.NorTA, ctrl.count.NorTA.spec, ctrl.compo.NorTA.spec,
     file = "../data/Data_by_NorTA_zeller.Rdata")

taxa.case.NorTA = colnames(case.count.NorTA); taxa.ctrl.NorTA = colnames(ctrl.count.NorTA)
save(taxa.case.NorTA, taxa.ctrl.NorTA, file = "../data/taxa_name_NorTA_zeller.Rdata")



