# phylogenetic hotspots


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
#library(arrow) # fast parquet file reading
library(data.table) # fast csv reading
library(castor) # fast tree reading
library(phyloregion) # PD calculations
library(raster)
library(sf)
library(ggplot2)
theme_set(theme_bw())
library(cowplot)
#library(rgdal)
library(beepr)



# Get phylogeny and distribution  ------------------------------------------

dist <- fread("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt")

# Remove non-accepted taxa from dist
nam <- fread("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", quote="")
table(dist$plant_name_id %in% nam$plant_name_id)
table(nam$taxon_status[nam$plant_name_id %in% dist$plant_name_id])
# no action required, unplaced will be sorted out by phylogeny matching for PD


# names(dist)
# table(dist$plant_name_id %in% t$tip.label)
# table(t$tip.label %in% dist$plant_name_id)

# remove doubtful location, extinct/introduced, and not included in phylogeny
# dist <- dist[dist$plant_name_id %in% t$tip.label,]
# rm(t)
dist <- dist[dist$extinct==0,]
dist <- dist[dist$introduced==0,]
dist <- dist[dist$location_doubtful==0,]
dist <- dist[!dist$area_code_l3=="",]
length(unique(dist$plant_name_id))

dist$area_code_l3 <- toupper(dist$area_code_l3)
length(unique(dist$area_code_l3))

# # get species for each bot country
# dist.list <- tapply(dist$plant_name_id, dist$area_code_l3, print)
# lengths(dist.list) # Species richness


# Get community matrix -----------------------------------------------------

#phyloregion::coldspots() / hotspots()
#phyloregion::long2dense()
#phyloregion::optimal_phyloregion()
#phyloregion::PD()
#phyloregion::PD_ses() # Phylogenetic diversity standardized for species richness
#phyloregion::phylo_endemism()

dist.mat <- long2sparse(dist, grids = "area_code_l3", species = "plant_name_id")
rm(dist)


# Get family count per country --------------------------------------------

# extract species IDs per bot country
res <- apply(dist.mat, 1, function(x){names(which(x == "1"))})
# get family names
res2 <- lapply(res, function(x){unique(nam$family[which(nam$accepted_plant_name_id %in% x)])})
#tmp <- data.frame(LEVEL3_COD = names(res2), family = lengths(res2)) # include list in df
saveRDS(res2, "families_per_level.rds")
rm(res,res2)

# Get phylogenies --------------------------------------------------------

# read in TACTed trees one at a time
phylonames <- dir("../DATA/phylos/TACT/", full.names = T)
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




# Null model scripts setup ---------------------------------------------

## outsource to run in parallel on the cluster for time efficiency ###
## write bash job array script and save data

### save data
# transfer to cluster manually, file too big for github

#check if all submat really the same:
#cn <- lapply(submat, colnames)
#for(i in 1:99)print(all(cn[i]%in%cn[i+1])) # all true, save space:
submat <- submat[[1]]

#save(list = c("submat", "subphy"), file="PD_nullmodel/comm_and_phy.RData")
# gdat = list(submat, subphy)
# saveRDS(gdat, file="PD_nullmodel/comm_and_phy.rds")


### build R scripts
Rscript <- "
  # # get rep number (=tree number)
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

print('PD_ses done')

# AvTD
tmp <- sparse2dense(submat[[rep]])
tmp <- as.data.frame(tmp)
rm(submat)
sphy <- subphy[[rep]]
rm(subphy)

res.avTD <- ecospat.calculate.pd(sphy, tmp, type='AvTD'', method='pairwise')
saveRDS(res.avTD, file=paste0('AvTD_, rep, '.rds'))

print('AvTD done'')

# TTD
res.TTD <- ecospat.calculate.pd(sphy, tmp, type='TTD', method='pairwise')
saveRDS(res.TTD, file=paste0('TTD_', rep, '.rds'))

print('TTD done')

all <- list(PD_ses_tipshuffle, res.avTD, res.TTD) #
saveRDS(all, file=paste0('indices_', rep, '.rds'))

ind$AvTD <- res.avTD
ind$TTD <- res.TTD
saveRDS(ind, file=paste0('ind_df_', rep, '.rds'))
"
cat(Rscript, file="PD_nullmodel/null_shuffle.R")

### build subordination bash script
bashscript2 <- '#!/bin/bash

#SBATCH --account PDiv
#SBATCH --job-name=nullshuffle
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=64gb
#SBATCH --cpus-per-task 1
#SBATCH --time 04:00:00
#SBATCH --output=nullshuffle.out

source ~/miniconda3/bin/activate R-env-4
Rscript null_shuffle.R $this_rep $this_model >logfile.txt
'
cat(bashscript2, file="PD_nullmodel/null_shuffle.sh")

### build bash array script
# The null model for separating patterns from processes and for
# contrasting against alternative hypotheses. Available null models include:
# “tipshuffle”: shuffles tip labels multiple times.
# “rowwise”: shuffles sites (i.e., varying richness) and keeping species
# occurrence frequency constant.
# “colwise”: shuffles species occurrence frequency and keeping site richness
# constant.
bashscript <- "#!/bin/bash
# submit_array.sh

#SBATCH --account PDiv
#SBATCH --job-name=nullshuffle_boss
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=1gb
#SBATCH --cpus-per-task 1
#SBATCH --time 00:05:00

rep=(1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72  73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100)
model=('tipshuffle' 'rowwise' 'colwise')


# pass on variables to child script
for ((m=0; m<=2; m++)) do
  this_model=${model[${m}]}
  export this_model
  for ((i=0; i<=99; i++)) do
    this_rep=${rep[${i}]}
    export this_rep
    sbatch null_shuffle.sh
  done
done
"
cat(bashscript, file="PD_nullmodel/null_shuffle_job_array.sh")

system({
"git add .
git commit -m 'nullmodel script changes'
git push"
})

# stuff on GenomeDK: git pull etc
# system2({
#   "ssh -T 'mtietje@login.genome.au.dk'"
# })
# #system("ssh -T 'mtietje@login.genome.au.dk'")


# Run my own nullmodel ----------------------------------------------------
 
# load("PD_nullmodel/comm_and_phy.RData")
# tree <- subphy[[1]]
# mat <- submat
# rm(subphy, submat)
# 
# # shuffle phylogeny tiplabels (ca 10 min for each iteration)
# set.seed(939)
# res <- matrix(ncol=100, nrow=nrow(mat))
# for(i in 1:100){
#   tree$tip.label <- sample(tree$tip.label)
#   res[,i] <- PD(mat, tree)
#   if(!i%%1)cat(i,"\r")
# }
# beep(2)
# # add observed
# #saveRDS(res, "res.rds")
# res <- cbind(res, PD(mat, subphy[[2]]))
# #PD_manual <- rowMeans(res)
# PD_sd_manual <- apply(res, 1, sd)
# 
# # plot histograms for each distribution and add observed PD
# resdf <- as.data.frame(res)
# names(resdf) <- c(paste0("rand_", 1:100), "PD_obs")
# saveRDS(resdf, "randomizedPD_100_and_obs.rds")
# 
# rest <- t(resdf)
# rest <- rest[rownames(rest)!="PD_obs",]
# resdf_long <- reshape2::melt(rest)
# bot <- unique(resdf_long$Var2)
# PD_obs <- PD(mat, subphy[[2]])
# PD_obs <- data.frame(PD_obs, Var2=names(PD_obs))
# tmp <- merge(resdf_long, PD_obs, all.x=TRUE)
# # add species richness
# s <- readRDS("fin_shp.rds")
# s <- st_drop_geometry(s)
# tmp <- merge(tmp, s[,c("LEVEL3_COD", "richness")], by.x="Var2", by.y="LEVEL3_COD", all.x=TRUE)
# tmp <- tmp[order(tmp$richness),]
# 
# ggplot(tmp[tmp$Var2 %in% bot[1:100],], aes(x=value, label=richness))+
#   geom_histogram(bins=120)+
#   facet_wrap(~Var2, ncol=5)+
#   geom_vline(aes(xintercept=PD_obs), col="red")+
#   geom_text(aes(x=0, y=75))
# 
# s <- s[!s$LEVEL3_COD=="BOU",]
# plot(tapply(tmp$value, tmp$Var2, sd), unique(tmp$PD_obs))
# plot(s$richness, res[,101]-rowMeans(res[,-101]), ylab="PDobs-PDrandom")
# plot(s$richness, (res[,101]-rowMeans(res[,-101]))/s$richness, ylab="(PDobs-PDrandom)/SR")

# sd does not increase the same as PD_obs. Small changes in small PD_obs are accompanied by big increases in SD, making it less likely the PDobs is outside the 2xSD - is this only a significance problem or also a ... z-score is the distance

# Repeated subsampling with n=minimum species number
## min number = 45 from the no.na dataset
## ____ Run this on cluster using bootstrap_job_array.sh ____________###
#submat <- submat[[1]]
#subphy <- subphy[[1]]
## all submats should be the same (given all tacted trees hold the same species)
# library(doMC)
# library(phytools)
# registerDoMC(cores = 2)
# set.seed(939)
# res <- matrix(ncol=100, nrow=nrow(submat)) # one column for each phylogeny
# res.sd <- matrix(ncol=100, nrow=nrow(submat)) # one column for each phylogeny
# # one phylo first, later all
# Sys.time()
# for(j in 1:100){
#   subphy <- subphy[[j]]
#   print(paste("Currently running phylogeny", j))
#   for(i in 1:nrow(submat)){ #
#     species.level <- colnames(submat)[which(submat[i,]==1)]
#     physub <- keep.tip(subphy, species.level)
#     tmp <- foreach(icount(1000)) %dopar% {
#       if(length(species.level)<45){
#         tmp <- NA
#         tmp
#       }else{
#         specs <- sample(species.level, 45)
#         physub2 <- keep.tip(physub, specs)
#         tmp <- sum(physub2$edge.length) # PD sensu stricto
#         tmp
#       }
#     }
#     res[i,j] <- mean(unlist(tmp), na.rm=T)
#     res.sd[i,j] <- sd(unlist(tmp), na.rm=T)
#     if(!i%%1)cat(i,"\r")
#   }-
# }
# Sys.time()
# beep(2)
# fil <- dir("PD_nullmodel", pattern="bootstrap_[1-9]", full.names = T)
# res <- lapply(fil, readRDS)
# 
# tmpPD <- lapply(res, FUN=function(x){x[[1]]})
# sdPD <- lapply(res, FUN=function(x){x[[2]]})
# pd.df <- as.data.frame(tmpPD)
# pdsd.df <- as.data.frame(sdPD)
# bootstrap.PD <- rowMeans(pd.df)
# bootstrap.PD.sd <- rowMeans(pdsd.df)
# 
# shp <- readRDS("fin_shp.rds")
# # remove BOU that has no data
# shp <- shp[!shp$LEVEL3_COD=="BOU",]
# 
# shp$bootstrapPD <- bootstrap.PD
# plot(shp$SES.PD, shp$bootstrapPD)
# plot(shp$PD_obs, shp$bootstrapPD)
# 
# shp$PD.SR <- shp$PD_obs/shp$richness
# shp$PD.AREA <- shp$PD_obs/shp$area
# ggplot(shp, aes(x=richness, y=PD_obs, col=log(area)))+
#   geom_point()+
#   facet_wrap(~CONTINENT)
# ggplot(shp, aes(x=richness, y=area))+
#   geom_point()+
#   facet_wrap(~CONTINENT, scales="free_y")
# ggplot(shp, aes(x=richness, y=SES.PD))+
#   geom_point()+
#   facet_wrap(~CONTINENT)
# # its actually still connected - just the other way now...
# ggplot(shp, aes(x=richness, y=PD.SR))+
#   geom_point()+
#   facet_wrap(~CONTINENT)
# ggplot(shp, aes(x=richness, y=log(PD.AREA)))+
#   geom_point()+
#   facet_wrap(~CONTINENT)
# ggplot(shp, aes(x=PD.AREA, y=PD_obs))+
#   geom_point()

# PD_manual <- rowMeans(res)
# PD_sd_manual <- apply(res, 1, sd)


# MPD + NRI (SES.MPD) ------------------------------------------------------
# library(PhyloMeasures)
# 
# system.time(load("PD_nullmodel/comm_and_phy.RData")) 
# 
# mat <- sparse2dense(submat)
# rm(submat)
# system.time(mpd <- mpd.query(subphy[[1]], mat, standardize = FALSE)) # 9sec
# 
# mpd <- matrix(ncol=100,nrow=nrow(mat))
# for(i in 1:100){
#   mpd[,i] <- mpd.query(subphy[[i]], mat, standardize = FALSE)
#   cat(i,"\r")
# }
# saveRDS(mpd, "mpd.rds")

# NRI runs on cluster with ses.mpd.R script
# ses.mpd <- mpd.query(subphy[[1]], mat, standardize = TRUE, null.model="uniform", reps=1000)







# Rarefaction a la Brody Sandel 2018 --------------------------------------

# Battling the assumption that species of different SR should have same
# environmental filtering, and PDI....
# One possible approach to removing richness-dependence is to rarefy all samples
# in the community matrix to the number of species in the smallest sample, then
# compute PDI.

load("PD_nullmodel/comm_and_phy.RData")
tree <- subphy[[1]]

submat <- mat[-as.numeric(which(rowSums(mat)<30)),] 
# remove bot countries with SR <30:
min(rowSums(submat))
# 41
PD_ses()

# rarefaction
res <- list(length=500)
rarefaction=41
no_zero <- apply(submat, 1, function(x)return(which(x==1))) # returns cols == 1

Sys.time()
set.seed(939)
for(i in 1:500){
  submat2 <- submat # reset matrix for each run
  no_zero_sub <- lapply(no_zero, sample, size=rarefaction)
  for(j in 1:nrow(submat2)){
    # set all not sampled entries to ZERO
    submat2[j, -as.numeric(no_zero_sub[[j]])] <- 0
  }
  res[[i]] <- PD_ses(submat2, tree, reps=100, "tipshuffle")
  if(!i%%1)cat(i,"\r")
}
Sys.time()

# stopped at 10h (=162 reps) for now

res.comb <- do.call(rbind, res)
# get original richness back in
res.comb$richness_original <- as.numeric(rep(rowSums(submat), length(res)))

ggplot(res.comb)+
  #geom_point(aes(y=zscore, x=richness_original, group=grids), alpha=.05)+
  stat_summary(aes(y=zscore, x=richness_original, group=grids), fun = "median", geom = "point")+
  ylab("median zscore rarefaction (n=162, sample size=41)")

ggplot(res.comb)+
  stat_summary(aes(y=zscore, x=richness_original, group=grids), fun = "sd", geom = "point")+
  ylab("SD zscore rarefaction (n=162, sample size=41)")

saveRDS(res.comb, "brody_rarefaction.rds")



# Botta-Dukat normalisaztion #####
# account for the non-normal distribution of the null-distribution
library(phyloregion)
library(ecospat)
load("PD_nullmodel/comm_and_phy.RData")
tree <- subphy[[1]]
sesPD_norm <- PD_ses_normal(submat, tree, reps=100, "tipshuffle")
sesPD_norm

saveRDS(sesPD_norm, "sesPD_norm.rds")






# SES.PD, AvTD, TTD --------------------------------------------------

# # SES.PD, AvTD, TTD
# pdnames <- dir("PD_nullmodel", pattern="ind_df", full.names = T)
# ind.list <- lapply(pdnames, readRDS)
# 
# #alt.list <- readRDS("PD_nullmodel/indices_2.rds") # 1=AvDT, 2=TTD
# #tmp <- lapply(ind.list, function(x){cbind(x[[1]], x[[2]])})
# 
# # ttd <- sapply(ind.list, "[[", 3)
# # TTD <- rowMeans(ttd)
# 
# tipshuffle.df <- data.frame(LEVEL3_COD = ind.list[[1]]$grids,
#                     richness = apply(sapply(ind.list, "[[", "richness"), 1, mean),
#                     PD_obs = apply(sapply(ind.list, "[[", "PD_obs"), 1, mean),
#                     pd_rand_mean = apply(sapply(ind.list, "[[", "pd_rand_mean"), 1, mean),
#                     pd_rand_sd = apply(sapply(ind.list, "[[", "pd_rand_sd"), 1, mean),
#                     pd_obs_rank = apply(sapply(ind.list, "[[", "pd_obs_rank"), 1, mean),
#                     SES.PD = apply(sapply(ind.list, "[[", "zscore"), 1, mean),
#                     pd_obs_p = apply(sapply(ind.list, "[[", "pd_obs_p"), 1, mean),
#                     reps = apply(sapply(ind.list, "[[", "reps"), 1, mean),
#                     AvTD = apply(sapply(ind.list, "[[", "AvTD"), 1, mean),
#                     TTD = apply(sapply(ind.list, "[[", "TTD"), 1, mean)
# )
# 
# # Careful with SD: from the function "pd_rand_sd <- apply(X = y, MARGIN = 2, FUN = var, na.rm = TRUE)" this is variance! --> zscore is fine: zscore <- (PD_obs - pd_rand_mean)/sqrt(pd_rand_sd)
# 
# # Removed entries with richness < 30? Zscore is not to be trusted: Central Limit Theorem
# table(tipshuffle.df$richness<30)
# tipshuffle.df[tipshuffle.df$richness<30, -c(1:3,9:11)] <- NA
# 
# # saveRDS(tipshuffle.df, "tipshuffle.rds")

# PD rowwise -------------------------------------------------------------
# pdnames <- dir("PD_nullmodel", pattern="rowwise", full.names = T)
# ind.list <- lapply(pdnames, readRDS)
# 
# #alt.list <- readRDS("PD_nullmodel/indices_2.rds") # 1=AvDT, 2=TTD
# #tmp <- lapply(ind.list, function(x){cbind(x[[1]], x[[2]])})
# 
# # ttd <- sapply(ind.list, "[[", 3)
# # TTD <- rowMeans(ttd)
# 
# rowshuffle.df <- data.frame(LEVEL3_COD = ind.list[[1]]$grids,
#                             richness = apply(sapply(ind.list, "[[", "richness"), 1, mean),
#                             PD_obs = apply(sapply(ind.list, "[[", "PD_obs"), 1, mean),
#                             pd_rand_mean = apply(sapply(ind.list, "[[", "pd_rand_mean"), 1, mean),
#                             pd_rand_sd = apply(sapply(ind.list, "[[", "pd_rand_sd"), 1, mean),
#                             pd_obs_rank = apply(sapply(ind.list, "[[", "pd_obs_rank"), 1, mean),
#                             SES.PD = apply(sapply(ind.list, "[[", "zscore"), 1, mean),
#                             pd_obs_p = apply(sapply(ind.list, "[[", "pd_obs_p"), 1, mean),
#                             reps = apply(sapply(ind.list, "[[", "reps"), 1, mean)
# )
# 
# # Careful with SD: from the function "pd_rand_sd <- apply(X = y, MARGIN = 2, FUN = var, na.rm = TRUE)" this is variance! --> zscore is fine: zscore <- (PD_obs - pd_rand_mean)/sqrt(pd_rand_sd)
# 
# # Removed entries with richness < 30? Zscore is not to be trusted: Central Limit Theorem
# table(rowshuffle.df$richness<30)
# rowshuffle.df[rowshuffle.df$richness<30, -c(1:3,9:11)] <- NA
# 
# # saveRDS(rowshuffle.df, "rowshuffle.rds")


# *** PE ------------------------------------------------------------------

# takes some time (couple minutes)
# load("PD_nullmodel/comm_and_phy.RData")
# tmp <- lapply(subphy, phylo_endemism, x=submat)
# PE <- colMeans(do.call(rbind,tmp))
# saveRDS(PE, "PE.rds")

# load results from GDK run
penames <- dir("PD_nullmodel", pattern="PE_", full.names = T)
ind.list <- lapply(penames, readRDS)

# alt.list <- readRDS("PD_nullmodel/indices_2.rds") # 1=AvDT, 2=TTD
# tmp <- lapply(ind.list, function(x){cbind(x[[1]], x[[2]])})

# ttd <- sapply(ind.list, "[[", )
# TTD <- rowMeans(ttd)

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

penames <- dir("PE_nullmodel", pattern="weighted=F", full.names = T)
ind.list <- lapply(penames, readRDS)
pe.df2 <- data.frame(LEVEL3_COD = ind.list[[1]]$grids,
                    richness = apply(sapply(ind.list, "[[", "richness"), 1, mean),
                    PE_obs = apply(sapply(ind.list, "[[", "PE_obs"), 1, mean),
                    pe_rand_mean = apply(sapply(ind.list, "[[", "pe_rand_mean"), 1, mean),
                    pe_rand_sd = apply(sapply(ind.list, "[[", "pe_rand_sd"), 1, mean),
                    pe_obs_rank = apply(sapply(ind.list, "[[", "pe_obs_rank"), 1, mean),
                    SES.PE = apply(sapply(ind.list, "[[", "zscore"), 1, mean),
                    pe_obs_p = apply(sapply(ind.list, "[[", "pe_obs_p"), 1, mean),
                    reps = apply(sapply(ind.list, "[[", "reps"), 1, mean)
)

# compare strict + weighted
plot(pe.df$PE_obs, pe.df2$PE_obs, xlab="weighted endemism (Rosauer et al. 2009)", ylab="strict endemism (Faith et al. 2004)")
cor.test(pe.df$PE_obs, pe.df2$PE_obs, method="s")
plot(pe.df$SES.PE, pe.df2$SES.PE)
cor.test(pe.df$SES.PE, pe.df2$SES.PE, method="s")
hist(pe.df2$SES.PE)

# Careful with SD: from the function "pd_rand_sd <- apply(X = y, MARGIN = 2, FUN = var, na.rm = TRUE)" this is variance! --> zscore is fine: zscore <- (PD_obs - pd_rand_mean)/sqrt(pd_rand_sd)

# Removed entries with richness < 30? Zscore is not to be trusted: Central Limit Theorem
which(pe.df$richness<30)
pe.df[pe.df$richness<30, -c(1:3,9:11)] <- NA



# *** Weighted endemism -------------------------------------------------------
# species richness inversely weighted by species ranges

WE <- weighted_endemism(submat) # one is enough
#rm(submat, subphy)






# RPD + RPE ---------------------------------------------------------------
# Mishler 2014:
# We also calculated for each grid cell two derived metrics that are new to this
# study and have been added as extensions of the Biodiverse package: RPD and
# RPE. Both of these indices are ratios that compare the PD and PE observed on
# the actual tree in the numerator to that observed on a comparison tree in the
# denominator. To make them easily comparable between analyses, the trees in
# both the numerator and the denominator are scaled such that branch lengths are
# calculated as a fraction of the total tree length. The comparison tree retains
# the actual tree topology but makes all branches of equal length. Thus, RPD is
# PD measured on the actual tree divided by PD measured on the comparison tree,
# while RPE is PE measured on the actual tree divided by PE measured on the
# comparison tree. In combination with the randomization test (below), this lets
# us examine the extent to which differential branch lengths matter to the
# patterns of PD and PE observed, which is important to our goals.

# install.packages("remotes")
# remotes::install_github("joelnitta/canaper")
library(canaper)
#data(phylocom)

load("PD_nullmodel/comm_and_phy.RData")
tree <- subphy[[1]]
rm(subphy)

# make data digestible for canaper (doesnt like ^numbers or -)
comm <- as.matrix(submat)
colnames(comm) <- paste0("sp", colnames(comm))
colnames(comm) <- gsub("-", "_", colnames(comm))
tree$tip.label <- paste0("sp", tree$tip.label)
tree$tip.label <- gsub("-", "_", tree$tip.label)

set.seed(071421)
Sys.time()
rand_test_results <- cpr_rand_test(
  comm, tree,
  null_model = "swap", 
  n_iterations=100
)
Sys.time()













# Assemble df + shapefile ------------------------------------------------

#rowshuffle.df <- readRDS("rowshuffle.rds")
tipshuffle.df <- readRDS("tipshuffle.rds")
## quick test
#plot(rowshuffle.df$PD_obs, tipshuffle.df$PD_obs)
#plot(rowshuffle.df$SES.PD, tipshuffle.df$SES.PD)
#names(rowshuffle.df)[4:8] <- paste0(names(rowshuffle.df)[4:8], "_RW")

mpd <- readRDS("mpd.rds")
mpd <- rowMeans(mpd)
ses.mpd <- readRDS("PD_nullmodel/ses_mpd.rds")
ses.mpd <- rowMeans(ses.mpd)
#PE <- readRDS("PE.rds")

s <- st_read("../DATA/shapefile_bot_countries/level3_fixed.gpkg")
s$LEVEL3_COD <- s$LEVEL_3_CO
s <- s[!s$LEVEL3_COD=="BOU",]
s <- merge(s, tipshuffle.df, by="LEVEL3_COD", all.x=TRUE)
#s <- merge(s, rowshuffle.df[,c(1,4:8)], by="LEVEL3_COD", all.x=TRUE)

# penalize SES.PD for richness?
plot(s$richness, s$SES.PD)
plot(s$richness, s$SES.PD/sqrt(s$richness))


s$mpd <- mpd
s$ses.mpd <- ses.mpd
s$WE <- WE

s$PE_obs <- pe.df$PE_obs
s$SES.PE <- pe.df$SES.PE
s$pe_obs_p <- pe.df$pe_obs_p


# *** Save  --------------------------------------------------------------

saveRDS(s, "fin_shape.rds")




# Plots -------------------------------------------------------------------


plot_grid(nrow=2,
          ggplot(tipshuffle.df, aes(x=richness, y=PD_obs, col=pd_obs_p<0.05))+
            geom_point()+
            scale_x_continuous("", trans="log", 
                               breaks = c(10, 100, 1000, 10000), limits = c(30, 23000))+
            #  scale_y_continuous("sqrt (PD)", trans="sqrt", breaks=c(100, 1000, 10000, 100000))+
            scale_color_discrete("diff from null dist", na.translate=F)+
            geom_line(aes(y=pd_rand_mean), col="grey")+
            geom_ribbon(aes(ymin=pd_rand_mean-2*sqrt(pd_rand_sd), 
                            ymax=pd_rand_mean+2*sqrt(pd_rand_sd)), alpha=0.3, color=NA)+
            theme(legend.position = c(x=0.15, y=0.85))
          ,
          ggplot(tipshuffle.df, aes(x=richness, y=SES.PD, col=pd_obs_p<0.05))+
            geom_point()+
            scale_x_continuous("", trans="log", limits = c(30, 23000))+
            theme(legend.position = "none")
          ,
          ggplot(tipshuffle.df, aes(x=richness, y=mpd))+
            geom_point()+
            scale_x_continuous("Species richness", trans="log", limits = c(30, 23000))
          ,
          ggplot(tipshuffle.df, aes(x=richness, y=ses.mpd))+
            geom_point()+
            scale_x_continuous("Species richness", trans="log", limits = c(30, 23000))
)

plot_grid(
ggplot(tipshuffle.df, aes(x=richness, y=PD_obs/richness, label=LEVEL3_COD))+
  geom_point(shape=NA)+
  geom_text(size=3)
,
ggplot(tipshuffle.df, aes(x=PD_obs, y=SES.PD, label=LEVEL3_COD))+
  geom_point(shape=NA)+
  geom_text(size=3)
)

# DOES THE Zscore SCALE WITH SR? yes
# ratio zscore to species richness:
plot(tipshuffle.df$richness, tipshuffle.df$SES.PD/tipshuffle.df$richness)

# Rank comparison SES.PD vs PD
PD_ranks <- order(tipshuffle.df$PD_obs, decreasing = T)
SES.PD_ranks <- order(tipshuffle.df$SES.PD, decreasing = T)
# biggest rank difference?
tipshuffle.df$LEVEL3_COD[which.max(PD_ranks - SES.PD_ranks)] # Wisconsin
rank.df = data.frame(rank.diff = PD_ranks - SES.PD_ranks, level3=s@data$LEVEL_NAME)

# lock factor levels order
rank.df <- rank.df[order(rank.diff),]
rank.df$level3 <- factor(rank.df$level3, levels = rank.df$level3)
pos <- as.numeric(rank.df$rank.diff<1)
pos[pos==0] <- -1

ggplot(rank.df, aes(x=level3, y=rank.diff, label=level3))+
  geom_bar(stat="identity")+
  geom_text(aes(y = pos*15,  angle = 45), size=3, hjust=as.numeric(rank.df$rank.diff>1))+
  facet_wrap(~factor(rank.diff<0), ncol = 1, scales = "free_x")#+
#  theme(axis.text.x = element_text(angle = 45))

# top.level3 <- tipshuffle.df$LEVEL3_COD[order(tipshuffle.df$PD_obs, decreasing = T)]
# top.level3.ses <- tipshuffle.df$LEVEL3_COD[order(abs(tipshuffle.df$SES.PD), decreasing = T)]
# data.frame(PD_order=top.level3, SES.PD_order=top.level3.ses)


