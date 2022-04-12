
  # get rep number (=tree number)
  args <- commandArgs()
  print(args)
  rep <- as.numeric(args[6])
  
  library(phyloregion) # PD calculations etc
  load('comm_and_phy.RData')
  PD_ses_tipshuffle <- PD_ses(submat[[1]], subphy[[1]], model='tipshuffle', reps=100) 
  saveRDS(PD_ses_tipshuffle, file=paste('tipshuffle_', rep, '.rds'))
