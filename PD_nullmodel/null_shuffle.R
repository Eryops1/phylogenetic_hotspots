
# get rep number (=tree number)
args <- commandArgs()
print(args)
rep <- as.numeric(args[6])
model <- args[7]

library(phyloregion) # PD calculations etc
library(ecospat)
load('comm_and_phy.RData')

# PD_ses
ind <- PD_ses(submat[[rep]], subphy[[rep]], model=model, reps=1000)
#saveRDS(PD_ses_tipshuffle, file=paste0(model,'_', rep, '.rds'))

print("PD_ses done")

# AvTD
tmp <- sparse2dense(submat[[rep]])
tmp <- as.data.frame(tmp)
rm(submat)
sphy <- subphy[[rep]]
rm(subphy)

res.avTD <- ecospat.calculate.pd(sphy, tmp, type="AvTD", method="pairwise")
#saveRDS(res.avTD, file=paste0("AvTD_", rep, '.rds'))

print("AvTD done")

# TTD
res.TTD <- ecospat.calculate.pd(sphy, tmp, type="TTD", method="pairwise")
#saveRDS(res.TTD, file=paste0("TTD_", rep, '.rds'))

print("TTD done")

all <- list(ind, res.avTD, res.TTD) #
saveRDS(all, file=paste0("indices_", rep, '.rds'))

ind$AvTD <- res.avTD
ind$TTD <- res.TTD
saveRDS(ind, file=paste0("ind_df_", rep, '.rds'))
