# Get PD endemism

setwd("~/Work/Aarhus/Science/PDiv")
rm(list = setdiff(ls(), lsf.str()))  
# library(data.table) # fast csv reading
# library(castor) # fast tree reading
library(phyloregion) # PD calculations

load("PD_nullmodel/comm_and_phy.RData")
dim(submat) # 330527
#nam = nam[nam$plant_name_id %in% colnames(submat),]
#length(unique(nam$genus))


# Get PD endemism -------------------------------------------------------

PDE_list = list()
Sys.time()
for(i in 1:length(subphy)){
  PDE_list[[i]] <- phylo_endemism(submat, subphy[[i]], weighted=F)
  # Strict endemism equates to the total amount of branch length found only in
  # the sample/s and is described by Faith et al. (2004)
  gc()
  cat(i, "\r")
}
Sys.time()


#saveRDS(PDE_list, "PDE_list.rds")


library(data.table)
tmp = as.data.table(PDE_list)
row.names(tmp) = names(PDE_list[[1]])
rowMeans(tmp)                    
