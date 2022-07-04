library(PhyloMeasures)
library(phyloregion)

load("comm_and_phy_sub.RData")

mat <- sparse2dense(submat)

ses.mpd <- matrix(ncol=100,nrow=nrow(mat))
for(i in 1:100){
  ses.mpd[,i] <- mpd.query(subphy[[i]], mat, standardize = TRUE, null.model="uniform", reps=1000)
  cat(i,"\r")
}

saveRDS(ses.mpd, "ses_mpd.rds")
