# phylogenetic hotspots


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
library(data.table) # fast csv reading
library(castor) # fast tree reading
library(phyloregion) # PD calculations
library(raster)
library(sf)
library(ggplot2)
theme_set(theme_bw())
library(cowplot)
library(beepr)




# Get distribution (WCVP data)  ------------------------------------------

dist <- fread("data/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt")
# this file is currently not published hence not in the repository

# Remove non-accepted taxa from dist
nam <- fread("data//wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", quote="")
table(dist$plant_name_id %in% nam$plant_name_id)
table(nam$taxon_status[nam$plant_name_id %in% dist$plant_name_id])
# no action required, unplaced will be sorted out by phylogeny matching for PD
nrow(nam[taxon_status=="Accepted" & taxon_rank=="Species" & genus_hybrid=="" & species_hybrid==""])


# remove doubtful location, extinct/introduced, clean LEVEL3 names
dist <- dist[dist$extinct==0,]
dist <- dist[dist$introduced==0,]
dist <- dist[dist$location_doubtful==0,]
dist <- dist[!dist$area_code_l3=="",]

dist$area_code_l3 <- toupper(dist$area_code_l3)
length(unique(dist$area_code_l3))


# Build community matrix -----------------------------------------------------

dist.mat <- long2sparse(dist, grids = "area_code_l3", species = "plant_name_id")
rm(dist)


# Get phylogenies --------------------------------------------------------

phylonames <- dir("TACT/", full.names = T)
phylonames <- phylonames[-grep("\\.sh$", phylonames)]
myfun <- function(x){return(read_tree(file=x, interpret_quotes = T))}
phylist <- lapply(phylonames, myfun)
rm(phylonames)

# match phylo and comm matrix
phylist2 <- lapply(phylist, match_phylo_comm, comm=dist.mat)
rm(phylist)

# check number of columns (should be identical)
unique(unlist(lapply(sapply(phylist2, "[[", "comm"), ncol)))
length(sapply(phylist2, "[[", "phy")[,1]$tip.label)

#rm(dist.mat)

submat <- lapply(phylist2, "[[", "comm")
subphy <- lapply(phylist2, "[[", "phy")
rm(phylist2)


submat <- submat[[1]]

save(list = c("submat", "subphy"), file="PD_nullmodel/comm_and_phy.RData")




## OUTSOURCE START ### ------
# 
# 
# run in parallel on the cluster for efficiency
# 
# submit jobarray to cluster for PDstd / PEstd. Takes each 4 hours and 64GB RAM.
# > null_shuffle_job_array.sh > null_shuffle.sh > null_shuffle.R
# Output: 100 files named ind_sf_[i].rds
# 
# repeat for PE switching to null_shuffle_PE.R 
# 
# 
## OUTSOURCE END ### --------







# Read data --------------------------------------------------
rm(subphy)

# SES.PD
pdnames <- dir("PD_nullmodel", pattern="ind_df", full.names = T)
ind.list <- lapply(pdnames, readRDS)

tipshuffle.df <- data.frame(LEVEL3_COD = ind.list[[1]]$grids,
                    richness = apply(sapply(ind.list, "[[", "richness"), 1, mean),
                    PD_obs = apply(sapply(ind.list, "[[", "PD_obs"), 1, mean),
                    pd_rand_mean = apply(sapply(ind.list, "[[", "pd_rand_mean"), 1, mean),
                    pd_rand_sd = apply(sapply(ind.list, "[[", "pd_rand_sd"), 1, mean),
                    pd_obs_rank = apply(sapply(ind.list, "[[", "pd_obs_rank"), 1, mean),
                    SES.PD = apply(sapply(ind.list, "[[", "zscore"), 1, mean),
                    pd_obs_p = apply(sapply(ind.list, "[[", "pd_obs_p"), 1, mean),
                    reps = apply(sapply(ind.list, "[[", "reps"), 1, mean),
                    AvTD = apply(sapply(ind.list, "[[", "AvTD"), 1, mean),
                    TTD = apply(sapply(ind.list, "[[", "TTD"), 1, mean)
)

# Removed entries with richness < 30: sample size too small

table(tipshuffle.df$richness<30)
tipshuffle.df[tipshuffle.df$richness<30, -c(1:3,9:11)] <- NA

saveRDS(tipshuffle.df, "data/sesPD.rds")



# Read data PE ------------------------------------------------------------

penames <- dir("PD_nullmodel", pattern="PE_", full.names = T)
ind.list <- lapply(penames, readRDS)

pe.df <- data.frame(LEVEL3_COD = ind.list[[1]]$grids,
                    richness = apply(sapply(ind.list, "[[", "richness"), 1, mean),
                    PE_obs = apply(sapply(ind.list, "[[", "PE_obs"), 1, mean),
                    pe_rand_mean = apply(sapply(ind.list, "[[", "pe_rand_mean"), 1, mean),
                    pe_rand_sd = apply(sapply(ind.list, "[[", "pe_rand_sd"), 1, mean),
                    pe_obs_rank = apply(sapply(ind.list, "[[", "pe_obs_rank"), 1, mean),
                    SES.PE = apply(sapply(ind.list, "[[", "zscore"), 1, mean),
                    pe_obs_p = apply(sapply(ind.list, "[[", "pe_obs_p"), 1, mean),
                    reps = apply(sapply(ind.list, "[[", "reps"), 1, mean)
)


# Removed entries with richness < 30
which(pe.df$richness<30)
pe.df[pe.df$richness<30, -c(1:3,9:11)] <- NA

saveRDS(pe.df, "data/sesPE.rds")


# Get weighted endemism -------------------------------------------------------

WE <- weighted_endemism(submat)



# Assemble df + shapefile ------------------------------------------------

pd.df <- readRDS("data/sesPD.rds")
pe.df <- readRDS("data/sesPE.rds")

s <- st_read("data/shapefile_bot_countries/level3_fixed.gpkg")
s$LEVEL3_COD <- s$LEVEL_3_CO
s <- merge(s, pd.df, by="LEVEL3_COD", all.x=TRUE)

s$WE <- WE

s$PE_obs <- pe.df$PE_obs
s$SES.PE <- pe.df$SES.PE
s$pe_obs_p <- pe.df$pe_obs_p


saveRDS(s, "data/fin_shape.rds")





# Botta-Dukat normalisaztion #####
# account for the non-normal distribution of the null-distribution
# library(phyloregion)
# library(ecospat)
# load("PD_nullmodel/comm_and_phy.RData")
# tree <- subphy[[1]]
# sesPD_norm <- PD_ses_normal(submat, tree, reps=100, "tipshuffle")
# sesPD_norm
# 
# saveRDS(sesPD_norm, "data/sesPD_norm.rds")
# 

