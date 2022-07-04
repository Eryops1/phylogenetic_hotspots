# ses.pd picante identiswap null model test


library(picante)

load("tmp.RData")

pic <- ses.pd(submat, phy, null.model = "independentswap",
              runs = 10, iterations = 10, include.root=TRUE)

saveRDS(pic, "picante_trial_run.rds")


