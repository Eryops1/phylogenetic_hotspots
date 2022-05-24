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
library(rgdal)
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
save(list = c("submat", "subphy"), file="PD_nullmodel/comm_and_phy.RData")


### build R scripts
Rscript <- "
  # get rep number (=tree number)
  args <- commandArgs()
  print(args)
  rep <- as.numeric(args[6])
  model <- args[7]
  
  library(phyloregion) # PD calculations etc
  load('comm_and_phy.RData')
  PD_ses_tipshuffle <- PD_ses(submat[[rep]], subphy[[rep]], model=model, reps=100)
  #PE_tipshuffle <- phylo_endemism((submat[[rep]], subphy[[rep]])
  saveRDS(PD_ses_tipshuffle, file=paste0(model,'_', rep, '.rds'))
"
cat(Rscript, file="PD_nullmodel/null_shuffle.R")

### build subordination bash script
bashscript2 <- '#!/bin/bash

#SBATCH --account PDiv
#SBATCH --job-name=nullshuffle
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=8gb
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00
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
# 
 load("PD_nullmodel/comm_and_phy.RData")



# 
# # shuffle phylogeny tiplabels (10 min or so for each iteration)
# set.seed(939)
# res <- matrix(ncol=100, nrow=nrow(submat[[1]]))
# for(i in 1:100){
#   subphy[[i]]$tip.label <- sample(subphy[[i]]$tip.label)
#   res[,i] <- PD(submat[[i]], subphy[[i]])
#   if(!i%%1)cat(i,"\r")
# }
# beep(2)
# PD_manual <- rowMeans(res)
# PD_sd_manual <- apply(res, 1, sd)


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
fil <- dir("PD_nullmodel", pattern="bootstrap_[1-9]", full.names = T)
res <- lapply(fil, readRDS)

tmpPD <- lapply(res, FUN=function(x){x[[1]]})
sdPD <- lapply(res, FUN=function(x){x[[2]]})
pd.df <- as.data.frame(tmpPD)
pdsd.df <- as.data.frame(sdPD)
bootstrap.PD <- rowMeans(pd.df)
bootstrap.PD.sd <- rowMeans(pdsd.df)

shp <- readRDS("fin_shp.rds")
# remove BOU that has no data
shp <- shp[!shp$LEVEL3_COD=="BOU",]

shp$bootstrapPD <- bootstrap.PD
plot(shp$SES.PD, shp$bootstrapPD)
plot(shp$PD_obs, shp$bootstrapPD)

shp$PD.SR <- shp$PD_obs/shp$richness
shp$PD.AREA <- shp$PD_obs/shp$area
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


# Pairwise distances ------------------------------------------------------
## NOT READY YET?
fil <- dir("PD_nullmodel", pattern="pairwise_dist_[1-9]", full.names = T)
res <- lapply(fil, readRDS)

tmppairD <- lapply(res, FUN=function(x){x[[1]]})
# sdpairD <- lapply(res, FUN=function(x){x[[2]]})
pairD.df <- as.data.frame(tmppairD)
# pairDsd.df <- as.data.frame(sdpairD)
pairD <- na.omit(rowMeans(pairD.df))


shp <- readRDS("fin_shp.rds")
# remove BOU that has no data
shp <- shp[!shp$LEVEL3_COD=="BOU",]
shp$pairD <- pairD




# Tipshuffle nullmodel ----------------------------------------------------

### read cluster results and get means
pdnames <- dir("PD_nullmodel", pattern="tipshuffle.*?.rds", full.names = T)

# 100 or 1000?
pdnames <- pdnames[grep("1000", pdnames)]

PD.list <- lapply(pdnames, readRDS)
tipshuffle.df <- data.frame(LEVEL3_COD = PD.list[[1]]$grids,
                    richness = apply(sapply(PD.list, "[[", "richness"), 1, mean),
                    PD_obs = apply(sapply(PD.list, "[[", "PD_obs"), 1, mean),
                    pd_rand_mean = apply(sapply(PD.list, "[[", "pd_rand_mean"), 1, mean),
                    pd_rand_sd = apply(sapply(PD.list, "[[", "pd_rand_sd"), 1, mean),
                    pd_obs_rank = apply(sapply(PD.list, "[[", "pd_obs_rank"), 1, mean),
                    SES.PD = apply(sapply(PD.list, "[[", "zscore"), 1, mean),
                    pd_obs_p = apply(sapply(PD.list, "[[", "pd_obs_p"), 1, mean),
                    reps = apply(sapply(PD.list, "[[", "reps"), 1, mean)
)
# names(tipshuffle.df)[-grep("LEVEL3|richness|reps", names(tipshuffle.df))] <- 
#   paste0(names(tipshuffle.df)[-grep("LEVEL3|richness|reps", names(tipshuffle.df))], "")

# Careful with SD: from the function "pd_rand_sd <- apply(X = y, MARGIN = 2, FUN = var, na.rm = TRUE)" this is variance! --> zscore is fine: zscore <- (PD_obs - pd_rand_mean)/sqrt(pd_rand_sd)

# Removed entries with richness < 30? Zscore is not to be trusted: Central Limit Theorem
table(tipshuffle.df$richness<30)
tipshuffle.df[tipshuffle.df$richness<30, -c(1,2)] <- NA


ggplot(tipshuffle.df, aes(x=richness, y=PD_obs, col=pd_obs_p<0.05))+
  geom_point()+
  scale_x_continuous("sqrt (species richness)", trans="sqrt", 
                     breaks = c(10, 100, 1000, 10000), limits = c(30, 23000))+
  scale_y_continuous("sqrt (PD)", trans="sqrt", breaks=c(100, 1000, 10000, 100000))+
  scale_color_discrete("diff from null dist", na.translate=F)+
  geom_line(aes(y=pd_rand_mean), col="grey")+
  geom_ribbon(aes(ymin=pd_rand_mean-2*sqrt(pd_rand_sd), 
              ymax=pd_rand_mean+2*sqrt(pd_rand_sd)), alpha=0.3, color=NA)+
  theme(legend.position = c(x=0.15, y=0.85), )

ggplot(tipshuffle.df, aes(x=richness, y=SES.PD, col=pd_obs_p<0.05))+
  geom_point()+
  scale_x_continuous(trans="log", limits = c(30, 23000))
ggplot(tipshuffle.df, aes(x=richness, y=pairD, col=pd_obs_p<0.05))+
  geom_point()+
  scale_x_continuous(trans="log", limits = c(30, 23000))


# DOES THE Zscore SCALE WITH SR??
# Does the SD of distributions scale with SES.PD?
ggplot(tipshuffle.df, aes(x=richness, y=sqrt(pd_rand_sd), col=pd_obs_p<0.05))+
  geom_point()#+
#  scale_x_continuous(limits = c(30, 23000))
## SD rises like crazy at the start and slows down for higher SR. An increase
## from 15k to 20k in SR only means and increase in SDrand of 100, whereas the
## same increase from 30 to 5k means an SD increase of 900.
ggplot(tipshuffle.df, aes(x=PD_obs, y=sqrt(pd_rand_sd), col=pd_obs_p<0.05))+
  geom_point()
## PD_randSD scales with observed PD non linearly 
ggplot(tipshuffle.df, aes(x=richness, y=sqrt(pd_rand_sd)/PD_obs))+
  geom_point()
## The ratio of PD_randSD : PD_obs gets rapidly Bigger for higher SR: It
## starts with 1:10 for small values and grows to 1:40 for high species
## richness. What does that mean: This does not affect the random PD curve
## position, but only its width: means that the random PD curve gets relatively
## narrower to the absolute value, so they can lay outside of this narrower
## distribution more easily?
# standard error of the mean equals the standard deviation divided by the square root of the sample size
ggplot(tipshuffle.df, aes(x=richness, y=sqrt(pd_rand_sd)/sqrt(richness)))+
  geom_point()
## Standard Error decreases more smoothly than SD

# Does scaling z score with SR again make sense?
ggplot(tipshuffle.df, aes(x=richness, y=SES.PD/richness, col=pd_obs_p<0.05))+
  geom_point()+
  scale_x_continuous(trans="log", limits = c(30, 23000))

# using SE instead of SD
tipshuffle.df$z.se=(tipshuffle.df$PD_obs-tipshuffle.df$pd_rand_mean)/(sqrt(tipshuffle.df$pd_rand_sd)/sqrt(tipshuffle.df$richness))
ggplot(tipshuffle.df, aes(x=PD_obs, y=z.se, col=pd_obs_p<0.05))+
  geom_point()
ggplot(tipshuffle.df, aes(x=PD_obs, y=SES.PD, col=pd_obs_p<0.05))+
  geom_point()


# # Rowwise nullmodel -------------------------------------------------------
# ### read cluster results and get means
# pdnames <- dir("PD_nullmodel", pattern="rowwise.*?.rds", full.names = T)
# rowwise.df <- data.frame(LEVEL3_COD = PD.list[[1]]$grids,
#                             richness = apply(sapply(PD.list, "[[", "richness"), 1, mean),
#                             PD_obs = apply(sapply(PD.list, "[[", "PD_obs"), 1, mean),
#                             pd_rand_mean = apply(sapply(PD.list, "[[", "pd_rand_mean"), 1, mean),
#                             pd_rand_sd = apply(sapply(PD.list, "[[", "pd_rand_sd"), 1, mean),
#                             pd_obs_rank = apply(sapply(PD.list, "[[", "pd_obs_rank"), 1, mean),
#                             SES.PD = apply(sapply(PD.list, "[[", "zscore"), 1, mean),
#                             pd_obs_p = apply(sapply(PD.list, "[[", "pd_obs_p"), 1, mean),
#                             reps = apply(sapply(PD.list, "[[", "reps"), 1, mean)
# )
# names(rowwise.df)[-grep("LEVEL3|richness|reps", names(rowwise.df))] <- 
#   paste0(names(rowwise.df)[-grep("LEVEL3|richness|reps", names(rowwise.df))], "_rw")
# 
# 
# 
# 
# # Colwise nullmodel -------------------------------------------------------
# ### read cluster results and get means
# pdnames <- dir("PD_nullmodel", pattern="colwise.*?.rds", full.names = T)
# PD.list <- lapply(pdnames, readRDS)
# colwise.df <- data.frame(LEVEL3_COD = PD.list[[1]]$grids,
#                          richness = apply(sapply(PD.list, "[[", "richness"), 1, mean),
#                          PD_obs = apply(sapply(PD.list, "[[", "PD_obs"), 1, mean),
#                          pd_rand_mean = apply(sapply(PD.list, "[[", "pd_rand_mean"), 1, mean),
#                          pd_rand_sd = apply(sapply(PD.list, "[[", "pd_rand_sd"), 1, mean),
#                          pd_obs_rank = apply(sapply(PD.list, "[[", "pd_obs_rank"), 1, mean),
#                          SES.PD = apply(sapply(PD.list, "[[", "zscore"), 1, mean),
#                          pd_obs_p = apply(sapply(PD.list, "[[", "pd_obs_p"), 1, mean),
#                          reps = apply(sapply(PD.list, "[[", "reps"), 1, mean)
# )
# names(colwise.df)[-grep("LEVEL3|richness|reps", names(colwise.df))] <- 
#   paste0(names(colwise.df)[-grep("LEVEL3|richness|reps", names(colwise.df))], "_cw")
# 
# 




# Assemble df + shapefile ------------------------------------------------

#dat <- merge.data.table(tipshuffle.df, rowwise.df, by=c("LEVEL3_COD", "richness", "reps"))
#dat <- merge.data.table(dat, colwise.df, by=c("LEVEL3_COD", "richness", "reps"))
#dat$PD_manual <- PD_manual
#dat$PD_sd_manual <- PD_sd_manual

#plot(dat$pd_rand_mean, dat$PD_manual)
#cor(dat$pd_rand_mean, dat$PD_manual)
#hist(dat$pd_rand_mean- dat$PD_manual)
# the phyloregion tipshuffle model is legit

s <- raster::shapefile("../DATA/shapefile_bot_countries/level3.shp")
s$LEVEL3_COD <- s$LEVEL_3_CO
s@data <- merge(s@data, tipshuffle.df, by="LEVEL3_COD", all.x=TRUE)


# *** PE ------------------------------------------------------------------

# takes some time (couple minutes)
load("PD_nullmodel/comm_and_phy.RData")
tmp <- mapply(phylo_endemism, submat, subphy)
tmp <- as.matrix(tmp)

pe <- data.frame(PE = apply(tmp, 1, mean), PE_sd = apply(tmp, 1, sd))
pe$LEVEL3_COD <- row.names(pe)
s@data <- merge(s@data, pe, all.x=TRUE)


# *** Weighted endemism -------------------------------------------------------
# species richness inversely weighted by species ranges

tmp <- mapply(weighted_endemism, submat[1]) # one is enough: distribution never changes
rm(submat, subphy)
tmp <- as.matrix(tmp)

we <- data.frame(WE = tmp[,1], LEVEL3_COD = row.names(tmp))
s@data <- merge(s@data, we, all.x=TRUE)



# *** Save  --------------------------------------------------------------


saveRDS(s, "fin_shape.rds")




