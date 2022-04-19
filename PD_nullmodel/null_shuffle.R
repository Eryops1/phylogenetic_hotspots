
  # get rep number (=tree number)
  args <- commandArgs()
  print(args)
  rep <- as.numeric(args[6])
  model <- args[7]
  
  library(phyloregion) # PD calculations etc
  load('comm_and_phy.RData')
  PD_ses_tipshuffle <- PD_ses(submat[[rep]], subphy[[rep]], model=model, reps=100) 
  saveRDS(PD_ses_tipshuffle, file=paste0(model,'_', rep, '.rds'))
