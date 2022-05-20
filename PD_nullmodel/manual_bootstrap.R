
# get rep number (=tree number)
args <- commandArgs()
print(args)
rep <- as.numeric(args[6])


library(phyloregion) # PD calculations etc
library(phytools)
library(castor)
load('comm_and_phy.RData')
submat <- submat[[rep]]
subphy <- subphy[[rep]]

res <- matrix(ncol=1, nrow=nrow(submat)) 
res.sd <- matrix(ncol=1, nrow=nrow(submat)) 
for(i in 1:nrow(submat)){ #
  species.level <- colnames(submat)[which(submat[i,]==1)]
  physub <- keep.tip(subphy, species.level)
  tmp <- for(j in 1:1000){ 
    if(length(species.level)<45){
      tmp <- NA
      tmp
    }else{
      specs <- sample(species.level, 45)
      physub2 <- get_subtree_with_tips(physub, specs) # is this the way to go? what about all other 
      tmp <- sum(physub2$edge.length) # PD sensu stricto
      tmp
      j <- j+1
    }
  }
  res[i,1] <- mean(unlist(tmp), na.rm=T)
  res.sd[i,1] <- sd(unlist(tmp), na.rm=T)
  if(!i%%1)cat(i,"\r")
}
  
res.all <- list(res, res.sd)  
saveRDS(res.all, file=paste0('bootstrap_', rep, '.rds'))
