
# get rep number (=tree number)
args <- commandArgs()
print(args)
rep <- as.numeric(args[6])
model <- args[7]

library(phyloregion) # PD calculations etc
library(ecospat)
load('comm_and_phy.RData')

# PD_ses
PD_ses_tipshuffle <- PD_ses(submat[[rep]], subphy[[rep]], model=model, reps=1000)
#PE_tipshuffle <- phylo_endemism((submat[[rep]], subphy[[rep]])
#saveRDS(PD_ses_tipshuffle, file=paste0(model,'_', rep, '.rds'))

print("PD_ses done")

# AvTD
tmp <- sparse2dense(submat[[rep]])
tmp <- as.data.frame(tmp)
res.avTD <- ecospat.calculate.pd(subphy[[rep]], tmp, type="AvTD")
#saveRDS(res.avTD, file=paste0("AvTD_", rep, '.rds'))

print("AvTD done")

# TTD
res.TTD <- ecospat.calculate.pd(subphy[[rep]], tmp, type="TTD")
#saveRDS(res.TTD, file=paste0("TTD_", rep, '.rds'))

print("TTD done")

all <- list(PD_ses_tipshuffle, res.avTD, res.TTD)
saveRDS(all, file=paste0("indices_", rep, '.rds'))




sp <- sample(subphy[[1]]$tip.label, 40000)
sphy <- keep.tip(subphy[[1]], sp)
tmp <- sparse2dense(submat[[1]])
tmp <- as.data.frame(tmp)
tmp <- tmp[,names(tmp) %in% sp]

system.time({
  pd <- ecospat.calculate.pd(sphy, tmp, method = "TTD")
})

