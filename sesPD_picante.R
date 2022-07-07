# ses.pd picante identiswap null model test


library(picante)

load("tmp.RData")
submat <- as.matrix(submat)


pic <- ses.pd(submat, phy, null.model = "independentswap",
              runs = 100, iterations = 100, include.root=TRUE)

saveRDS(pic, "picante_run100.rds")


## results
# 
# dat <- readRDS("picante_trial_run.rds")
