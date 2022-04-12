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



# dat <- fread("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", quote = "")
# dist <- fread("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt")
# write_parquet(dat, "../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_names.parquet")
# write_parquet(dist, "../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.parquet")

# dat <- read_parquet("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_names.parquet")
# dist <- read_parquet("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.parquet")



# Get phylogeny and distribution  ------------------------------------------

t <- read_tree(file="../DATA/phylos/seed_plant_TOL_WCVP_IDs.tre")
dist <- fread("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt")

names(dist)
table(dist$plant_name_id %in% t$tip.label)
table(t$tip.label %in% dist$plant_name_id)

# remove doubtful location, extinct/introduced, and not included in phylogeny
dist <- dist[dist$plant_name_id %in% t$tip.label,]
rm(t)
dist <- dist[dist$extinct==0,]
dist <- dist[dist$introduced==0,]
dist <- dist[dist$location_doubtful==0,]
length(unique(dist$plant_name_id))

dist$area_code_l3 <- toupper(dist$area_code_l3)
length(unique(dist$area_code_l3))

# get species for each bot country
#dist.list <- tapply(dist$plant_name_id, dist$area_code_l3, print)
#lengths(dist.list) # Species richness


# Get community matrix -----------------------------------------------------

#phyloregion::coldspots() / hotspots()
phyloregion::long2dense()
phyloregion::optimal_phyloregion()
phyloregion::PD()
phyloregion::PD_ses()
phyloregion::phylo_endemism()

dist.mat <- long2sparse(dist, grids = "area_code_l3", species = "plant_name_id")
rm(dist)

# Get PD etc ---------------------------------------------------------------

t <- read_tree(file="../DATA/phylos/seed_plant_TOL_WCVP_IDs.tre")
s <- raster::shapefile("../DATA/wgsrpd-master/level3/level3.shp")

# match phylo and comm matrix
subphy <- match_phylo_comm(t, dist.mat)$phy
submat <- match_phylo_comm(t, dist.mat)$com

length(colnames(submat))
length(subphy$tip.label)

rm(t, dist.mat)

s@data$PD <- PD(submat, subphy)
#s@data$PD_ses <- PD_ses(submat, subphy)
s@data$PE <- phylo_endemism(submat, subphy)
s@data$SR <- as.numeric(rowSums(submat))

PD_ses(submat, subphy)

# format to sf object for plotting
shp <- st_as_sf(s)



# Maps --------------------------------------------------------------------

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

ggplot(shp, aes(PE, PD))+
  geom_point()+
  scale_y_sqrt()+
  scale_x_sqrt()



