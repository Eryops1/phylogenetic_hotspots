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

# match phylso and comm matrix
phylist2 <- lapply(phylist, match_phylo_comm, comm=dist.mat)
rm(phylist)

# check number of columns (should be identical)
unique(unlist(lapply(sapply(phylist2, "[[", "comm"), ncol)))
length(sapply(phylist2, "[[", "phy")[,1]$tip.label)

rm(dist.mat)


# Assemble shapefile ------------------------------------------------------

s <- raster::shapefile("../DATA/wgsrpd-master/level3/level3.shp")



# *** PD ------------------------------------------------------------------

submat <- lapply(phylist2, "[[", "comm")
subphy <- lapply(phylist2, "[[", "phy")
rm(phylist2)

# # takes some time (couple minutes) - can also delete this one and use null model setting for observed PD
# tmp <- mapply(PD, submat, subphy)
# tmp <- as.matrix(tmp)
# 
# s@data$PD <- apply(tmp, 1, mean)
# s@data$PD_sd <- apply(tmp, 1, sd)
# #s@data$PD_sd <- apply(tmp, 1, se)



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



# *** Null model PD -------------------------------------------------------

# PD_ses --> Phylogenetic diversity standardized for species richness

## outsource to run in parallel on the cluster for time efficiency ###
## write bash job array script and save data

### save data
save(list = c("submat", "subphy"), file="comm_and_phy.RData")

### build R script
Rscript <- "
  # get rep number (=tree number)
  args <- commandArgs()
  print(args)
  rep <- as.numeric(args[6])
  
  library(phyloregion) # PD calculations etc
  load('comm_and_phy.RData')
  PD_ses_tipshuffle <- PD_ses(submat[[rep]], subphy[[rep]], model='tipshuffle', reps=100) 
  saveRDS(PD_ses_tipshuffle, file=paste0('tipshuffle_', rep, '.rds'))
"
cat(Rscript, file="tipshuffle.R")

### build bash array script
bashscript <- "#!/bin/bash
# submit_array.sh

#SBATCH --account PDiv
#SBATCH --job-name=tipshuffle_master
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=1gb
#SBATCH --cpus-per-task 1
#SBATCH --time 00:05:00

rep=seq 1 67

export rep
# write bash scripts
for ((i=0; i<=66; i++)) do
  this_rep=${rep[${i}]}
  export this_rep
  sbatch tipshuffle.sh
done
"
cat(bashscript, file="tipshuffle_job_array.sh")

### build subordination bash script
bashscript2 <- '#!/bin/bash

#SBATCH --account PDiv
#SBATCH --job-name=tipshuffle
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=8gb
#SBATCH --cpus-per-task 1
#SBATCH --time 00:15:00
#SBATCH --output=/dev/null

source ~/miniconda3/bin/activate R-env-4
Rscript tipshuffle.R $this_rep > log_"$SLURM_JOB_ID"_"$this_rep".txt
'
cat(bashscript2, file="tipshuffle.sh")





### read cluster results and get means
pdnames <- dir("../DATA/PDiv/PD_nullmodel", full.names = T)
PD.list <- lapply(pdnames, readRDS)
PD.df <- data.frame(PD_obs_mean <- apply(sapply(PD.list, "[[", "PD_obs"), 1, mean),
)


s@data <- merge(s@data, PD_ses_tipshuffle, by="LEVEL3_COD")

PD_ses_rowwise <- PD_ses(submat[[1]], subphy[[1]], model="rowwise", reps=100)
PD_ses_colwise <- PD_ses(submat[[1]], subphy[[1]], model="colwise", reps=100)

plot(PD_ses_tipshuffle)

# *** Save shapefile ------------------------------------------------------

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
