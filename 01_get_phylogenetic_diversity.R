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

# match phylso and comm matrix
phylist2 <- lapply(phylist, match_phylo_comm, comm=dist.mat)
rm(phylist)

# check number of columns (should be identical)
unique(unlist(lapply(sapply(phylist2, "[[", "comm"), ncol)))
length(sapply(phylist2, "[[", "phy")[,1]$tip.label)

rm(dist.mat)

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
#SBATCH --output=/dev/null

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
for ((m=0; m<=2; i++)) do
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

# git add commit push --> GenomeDK


# Tipshuffle nullmodel ----------------------------------------------------

### read cluster results and get means
pdnames <- dir("PD_nullmodel", pattern="tipshuffle.*?.rds", full.names = T)
PD.list <- lapply(pdnames, readRDS)
PD.df <- data.frame(LEVEL3_COD = PD.list[[1]]$grids,
                    richness = apply(sapply(PD.list, "[[", "richness"), 1, mean),
                    PD_obs_mean = apply(sapply(PD.list, "[[", "PD_obs"), 1, mean),
                    pd_rand_m_mean = apply(sapply(PD.list, "[[", "pd_rand_mean"), 1, mean),
                    pd_rand_sd_mean = apply(sapply(PD.list, "[[", "pd_rand_sd"), 1, mean),
                    pd_obs_rank_mean = apply(sapply(PD.list, "[[", "pd_obs_rank"), 1, mean),
                    zscore_mean = apply(sapply(PD.list, "[[", "zscore"), 1, mean),
                    pd_obs_p_mean = apply(sapply(PD.list, "[[", "pd_obs_p"), 1, mean),
                    reps = apply(sapply(PD.list, "[[", "reps"), 1, mean)
)
names()

s@data <- merge(s@data, PD.df, by="LEVEL3_COD", all.x=TRUE)

plot(s@data[,-grep("LEVEL|reps", names(s@data))])
hist(s@data$zscore_mean)

ggplot(data=s@data, aes(x=richness, y=PD_obs_mean))+
  geom_point()

ggplot(data=s@data, aes(x=pd_rand_m_mean, y=PD_obs_mean, col=log(richness)))+
  geom_point()+
  geom_abline(x=1)+
  scale_x_continuous(trans="sqrt")+
  scale_y_continuous(trans="sqrt")
# observed PD is usually lower than expected

ggplot(data=s@data, aes(x=richness, y=zscore_mean))+
  geom_point()+
  geom_hline(yintercept=-1)+
  scale_x_continuous(trans = "sqrt")
# zscore becomes more negative with increasing richness: less species than expected

ggplot(data=s@data, aes(x=richness, y=PD_obs_mean))+
  geom_point()+
  geom_point(aes(y=pd_rand_m_mean), col="salmon")+
  scale_x_continuous(trans="sqrt")+
  scale_y_continuous(trans="sqrt")

hist(s@data$PD_obs_mean/s@data$richness)


# PD_obs: observed PD in community
# pd_rand_mean: mean PD in null communities
# pd_obs_rank: Rank of observed PD vs. null communities
# pd_obs_z: Standardized effect size of PD vs. null communities = (PD_obs - pd_rand_mean) / pd_rand_sd
# pd_obs_p: P-value (quantile) of observed PD vs. null communities = mpd_obs_rank / iter + 1


# Rowwise nullmodel -------------------------------------------------------
### read cluster results and get means
pdnames <- dir("PD_nullmodel", pattern="rowwise.*?.rds", full.names = T)
PD.list <- lapply(pdnames, readRDS)
PD.df <- data.frame(LEVEL3_COD = PD.list[[1]]$grids,
                    richness = apply(sapply(PD.list, "[[", "richness"), 1, mean),
                    PD_obs_mean = apply(sapply(PD.list, "[[", "PD_obs"), 1, mean),
                    pd_rand_m_mean = apply(sapply(PD.list, "[[", "pd_rand_mean"), 1, mean),
                    pd_rand_sd_mean = apply(sapply(PD.list, "[[", "pd_rand_sd"), 1, mean),
                    pd_obs_rank_mean = apply(sapply(PD.list, "[[", "pd_obs_rank"), 1, mean),
                    zscore_mean = apply(sapply(PD.list, "[[", "zscore"), 1, mean),
                    pd_obs_p_mean = apply(sapply(PD.list, "[[", "pd_obs_p"), 1, mean),
                    reps = apply(sapply(PD.list, "[[", "reps"), 1, mean)
)

ggplot(data=PD.df, aes(x=richness, y=pd_rand_m_mean))+
  geom_point()
ggplot(data=s@data, aes(x=richness, y=pd_rand_m_mean))+
  geom_point()





# Colwise nullmodel -------------------------------------------------------
### read cluster results and get means
pdnames <- dir("PD_nullmodel", pattern="colwise.*?.rds", full.names = T)
PD.list <- lapply(pdnames, readRDS)
PD.df <- data.frame(LEVEL3_COD = PD.list[[1]]$grids,
                    richness = apply(sapply(PD.list, "[[", "richness"), 1, mean),
                    PD_obs_mean = apply(sapply(PD.list, "[[", "PD_obs"), 1, mean),
                    pd_rand_m_mean = apply(sapply(PD.list, "[[", "pd_rand_mean"), 1, mean),
                    pd_rand_sd_mean = apply(sapply(PD.list, "[[", "pd_rand_sd"), 1, mean),
                    pd_obs_rank_mean = apply(sapply(PD.list, "[[", "pd_obs_rank"), 1, mean),
                    zscore_mean = apply(sapply(PD.list, "[[", "zscore"), 1, mean),
                    pd_obs_p_mean = apply(sapply(PD.list, "[[", "pd_obs_p"), 1, mean),
                    reps = apply(sapply(PD.list, "[[", "reps"), 1, mean)
)

ggplot(data=PD.df, aes(x=richness, y=pd_rand_m_mean))+
  geom_point()
ggplot(data=s@data, aes(x=richness, y=pd_rand_m_mean))+
  geom_point()





# Assemble shapefile ------------------------------------------------------

# PD_ses --> Phylogenetic diversity standardized for species richness

s <- raster::shapefile("../DATA/wgsrpd-master/level3/level3.shp")

# *** PE ------------------------------------------------------------------

# takes some time (couple minutes)
tmp <- mapply(phylo_endemism, submat, subphy)
tmp <- as.matrix(tmp)

pe <- data.frame(PE = apply(tmp, 1, mean), PE_sd = apply(tmp, 1, sd))
pe$LEVEL3_COD <- row.names(pe)
s@data <- merge(s@data, pe, all.x=TRUE)


# *** SR ------------------------------------------------------------------
sr <- data.frame(SR=rowSums(submat[[1]]))
sr$LEVEL3_COD <- row.names(sr)
s@data <- merge(s@data, sr, all.x=TRUE)


# *** Save  --------------------------------------------------------------

s@data <- merge(s@data, PD.df, by="LEVEL3_COD", all.x=TRUE)

saveRDS(s, "fin_shape.rds")









# Maps --------------------------------------------------------------------

# format to sf object for plotting
shp <- st_as_sf(s)

# transform to Behrmann projection
shp <- st_transform(shp, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs") # Behrmann

(pd_map <- ggplot(shp) + 
   geom_sf(aes(fill=PD),lwd=0, col=NA) + 
   #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
   #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
  #                        begin = lcol, end = sqrt(ucol))+
   scale_fill_viridis_c("PD", option = "plasma")+ #, 
   theme(legend.position = c(0.18, 0.3),
         legend.key.height = unit(6,"mm"),
         legend.background = element_blank(),
         legend.key = element_blank(),
         panel.background = element_blank(),
         #panel.border = element_blank(),
         text = element_text(size = 10),
         # axis.ticks.x = element_line(color="white"),
         # axis.text.x = element_text(color="white"),
         # plot.margin = margin(0, 0, 0, 0.5, "cm")
   )+
   #coord_sf(expand = F, label_axes = "-N--", ylim=c(-6200000, 8200000))+
   xlab(" ")
)
(pd_sd_map <- ggplot(shp) + 
    geom_sf(aes(fill=PD_sd),lwd=0, col=NA) + 
    #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
    #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
    #                        begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("PD_sd", option = "plasma")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")
)
(pd_se_map <- ggplot(shp) + 
    geom_sf(aes(fill=PD_sd/sqrt(SR)),lwd=0, col=NA) + 
    #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
    #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
    #                        begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("PD_se", option = "plasma")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")
)

# (pd_ses_map <- ggplot(shp) + 
#     geom_sf(aes(fill=PD_ses),lwd=0, col=NA) + 
#     #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
#     #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
#     #                        begin = lcol, end = sqrt(ucol))+
#     scale_fill_viridis_c("PD_ses", option = "plasma")+ #, 
#     theme(legend.position = c(0.18, 0.3),
#           legend.key.height = unit(6,"mm"),
#           legend.background = element_blank(),
#           legend.key = element_blank(),
#           panel.background = element_blank(),
#           #panel.border = element_blank(),
#           text = element_text(size = 10),
#           # axis.ticks.x = element_line(color="white"),
#           # axis.text.x = element_text(color="white"),
#           # plot.margin = margin(0, 0, 0, 0.5, "cm")
#     )+
#     #coord_sf(expand = F, label_axes = "-N--", ylim=c(-6200000, 8200000))+
#     xlab(" ")
# )

(pe_map <- ggplot(shp) + 
    geom_sf(aes(fill=PE),lwd=0, col=NA) + 
    #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
    #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
    #                        begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("PE", option = "plasma")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          #panel.border = element_blank(),
          text = element_text(size = 10),
          # axis.ticks.x = element_line(color="white"),
          # axis.text.x = element_text(color="white"),
          # plot.margin = margin(0, 0, 0, 0.5, "cm")
    )+
    #coord_sf(expand = F, label_axes = "-N--", ylim=c(-6200000, 8200000))+
    xlab(" ")
)

(sr_map <- ggplot(shp) + 
    geom_sf(aes(fill=SR),lwd=0, col=NA) + 
    #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
    #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
    #                        begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("SR", option = "plasma")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          #panel.border = element_blank(),
          text = element_text(size = 10),
          # axis.ticks.x = element_line(color="white"),
          # axis.text.x = element_text(color="white"),
          # plot.margin = margin(0, 0, 0, 0.5, "cm")
    )+
    #coord_sf(expand = F, label_axes = "-N--", ylim=c(-6200000, 8200000))+
    xlab(" ")
)

plot_grid(pd_map, pe_map, sr_map, nrow = 2)

plot_grid(nrow=2,
ggplot(shp, aes(PE, PD))+
  geom_point()+
  scale_y_sqrt()+
  scale_x_sqrt()
,
ggplot(shp, aes(SR, PD))+
  geom_point()+
  scale_y_sqrt()+
  scale_x_sqrt()
,
ggplot(shp, aes(SR, PE))+
  geom_point()+
  scale_y_sqrt()+
  scale_x_sqrt()
)
