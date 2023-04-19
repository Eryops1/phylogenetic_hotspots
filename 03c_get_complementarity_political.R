# phylogenetic hotspots


rm(list = setdiff(ls(), lsf.str()))  
library(data.table) # fast csv reading
library(castor) # fast tree reading
library(phyloregion) # PD calculations
library(raster)
library(sf)
library(ggplot2)
theme_set(theme_classic())
library(cowplot)
library(terra)

# Load data

# only if still calculating things
load("PD_nullmodel/comm_and_phy.RData")
source("99_functions.R")


pc <- st_read("data/political_countries.gpkg")
shp <- st_read("data/shapefile_bot_countries/level3_fixed.gpkg")

# get the merged countries
#sf::st_intersection(pc[,1], shp[,1])
# countries = c("CAN", "USA", "MXC", "BRZ", "AGE", "CLC", "SAF", "RSA", "AUS", "CHN", "MOR")
# res = list()
# for(i in 1:length(countries)){
#   res[[i]] <- sf::st_intersection(pc[pc$LEVEL_3_CO==countries[i],1], shp[,1])
# }
# mc = sapply(res, "[[", "LEVEL_3_CO")
# names(mc) = countries
# cdt = unlist(mc)
# cdt = data.table(LEVEL3_COD_merge = names(cdt), LEVEL3_COD = cdt)
# cdt$LEVEL3_COD_merge = gsub("[0-9]{1,2}", "", cdt$LEVEL3_COD_merge)
# cdt
# 
# # merge the presence absence matrix
# pos = match(row.names(submat), cdt$LEVEL3_COD) # position of first argument in second, etc
# 
# row.names(submat)[which(!is.na(pos))] = cdt$LEVEL3_COD_merge[na.omit(pos)]
# submat
# 
# # merge rows with identical names
# mat = matrix(nrow=length(countries), ncol=ncol(submat))
# for(i in 1:length(countries)){
#   tmp = submat[row.names(submat) %in% countries[i],]
#   tn = colSums(tmp)
#   tn[tn>1] = 1
#   mat[i, ] = tn
# }
# 
# row.names(mat) = countries
# 
# submat2 = submat[!row.names(submat) %in% countries,]
# row.names(submat2)
# submat_store = submat
# ## actual merge
# submat = rbind(submat2, mat)
# row.names(submat)
# 
# saveRDS(submat, "data/merged_distribution_matrix.rds")
submat <- readRDS("data/merged_distribution_matrix.rds")

# Get SR complementarity -----------------------------------------------------

# greed <- sr_greedy(submat, m=0) #, n=10
# saveRDS(greed, "data/sr_complementarity_political.rds")

shp <- st_read("data/political_countries.gpkg")
shp$LEVEL3_COD = shp$LEVEL_3_CO

greed = readRDS("data/sr_complementarity_political.rds")
res <- data.table(LEVEL3_COD = names(greed[[1]]),
                  sr_complementarity = greed[[1]])

tmp = data.table(sr_added = greed[[2]],
                 LEVEL3_COD = greed[[3]])
res <- merge(res, tmp, all=T)
shp <- merge(shp, res, all.x=T)


# Get PD complementarity -----------------------------------------------------

# res <- list()
# for(i in 1:1){
#   greed <- pd_greedy(submat, phylo=subphy[[i]])
#   res[[i]] <- greed
#   cat(i, "\r")
# }
# #saveRDS(res, "data/pd_complementarity_political.rds")
# gc()

# read results:
fnames <- dir("data/political_countries/")[grep("pd_complementarity_political",
                                                dir("data/political_countries"))]
lpd <- list()
for(i in 1:100){
  lpd[[i]] <- readRDS(paste0("data/political_countries/", fnames[i]))
}

pd_complementarity = sapply(lpd, "[[", 1)
rowSums(pd_complementarity, na.rm=T)
# this is the same everywhere, across all TACT trees countries either contribute
# or not (not=36)
res <- data.table(pd_complementarity = pd_complementarity[,1],
                  LEVEL3_COD = row.names(pd_complementarity))

pd_added = sapply(lpd, "[[", 2)
pd_area = sapply(lpd, "[[", 3)

# Datatable of dimensions col=100 tress, row=country position. If there are more
# than 1 country per row, that means the position of this country is not always
# the same. Quantify:
apply(pd_area, 1, function(x)length(unique(x)))

# the order is not always the same due to minor differences between the trees.
# Get the most common country for each position that is selected by the
# algorithm:
#most_common <- apply(pd_area, 1, function(x)names(sort(table(x),decreasing=TRUE))[1])
# this does not work as i have double or even trippe countries and some even missing.

# get average value for each country instead, ignore order
tmp = data.table(pd_added = as.numeric(pd_added),
                 LEVEL3_COD = as.character(pd_area))
average_PD_added <- tapply(tmp$pd_added, tmp$LEVEL3_COD, mean)
res2 <- data.table(LEVEL3_COD = names(average_PD_added),
                   pd_added = average_PD_added)
res <- merge(res, res2, all.x=T, by="LEVEL3_COD")
shp <- merge(shp, res, all.x=T)


# Get PD endemism complementarity -------------------------------------------

# PDendemism is already its complementarity value, also the algirthm does not
# work since it alters absences of already covered countries, artificially
# increasing other countries endemism values. Use raw PDE values, ordered.

# Get PD endmism for political countries
pde <- readRDS("data/political_countries/PDE_list_political.rds")
tmp = as.data.table(pde)
row.names(tmp) = names(pde[[1]])
PDE = rowMeans(tmp)  
names(PDE) = row.names(tmp)
tmp <- data.table(PDE=PDE, LEVEL3_COD=names(PDE))


shp <- merge(shp, tmp, all.x=T)
shp$pde_complementarity <- shp$PDE > 0

#View(st_drop_geometry(shp[,c("LEVEL3_COD", "LEVEL_NAME", "PDE", "PE_obs",
#                             "pde_complementarity", "PD_obs")]))











# Get number of countries for wanted PERCENTAGES captured --------------------

# percentages
shp$sr_added_perc = shp$sr_added/sum(shp$sr_added, na.rm=T)
shp$pd_added_perc = shp$pd_added/sum(shp$pd_added, na.rm=T)
shp$pde_added_perc = shp$PDE/sum(shp$PDE, na.rm=T)

# set zeros in PDE to NA
shp$pde_added_perc[shp$pde_added_perc==0] = NA

perc_sr = c()
for(i in seq(.1, 1, by=.1)){
  shpsr = na.omit(shp$sr_added_perc)
  tmp = c()
  while(sum(tmp) <= i & length(shpsr)>0) {
    # take percentage and add to set
    tmp <- c(tmp, max(shpsr, na.rm=T))
    shpsr = shpsr[-which.max(shpsr)]
  }
  perc_sr = c(perc_sr, length(tmp))
  cat(i, "\r")
}
perc_sr
# To conserve 50% of SR, we need 9 countries

perc_pd = c()
for(i in seq(.1, 1, by=.1)){
  shppd = na.omit(shp$pd_added_perc)
  tmp = c()
  while(sum(tmp) <= i & length(shppd)>0) {
    # take percentage and add to set
    tmp <- c(tmp, max(shppd, na.rm=T))
    shppd = shppd[-which.max(shppd)]
  }
  perc_pd = c(perc_pd, length(tmp))
  cat(i, "\r")
}
perc_pd
# To conserve 50% of PD, we need 22 countries

perc_pde = c()
for(i in seq(.1, 1, by=.1)){
  shppde = shp$pde_added_perc
  tmp = c()
  while(sum(tmp) <= i & length(shppde)>0) {
    # take percentage and add to set
    tmp <- c(tmp, max(shppde, na.rm=T))
    shppde = shppde[-which.max(shppde)]
  }
  perc_pde = c(perc_pde, length(tmp))
  cat(i, "\r")
}
perc_pde
# To conserve 50% of PDE, we need 6 countries


# Barplot (Fig 4)
pdt = data.table(percent_captured = seq(.1, 1, by=.1),
           SR = perc_sr, 
           PD = perc_pd, 
           PD_endemism = perc_pde)
pdt = melt(pdt, id="percent_captured")

bc <- c("#35abc4", "#4b9e31", "#eeea40")

(barplot = ggplot(pdt, aes(x= percent_captured, y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_discrete("")+
  ylab("Political countries")+
  xlab("Proportion captured")+
  scale_fill_manual("", values=c("grey50", bc[1], bc[2]), labels=c('SR', 'PD', 'PD endemism'))+
  theme(legend.position=c(.25,.9), 
        legend.background=element_blank())
)
ggsave("figures/barplot_percentages_political_countries.png", width=3.5, height=3.5, dpi=600, bg="white")
  






## Define top 7 / 50% -------------------------------------------------------


# SR top 10 and 50%
shp <- shp[order(shp$sr_added, decreasing=T), ]
shp$sr_comp_top10 <- shp$sr_complementarity
shp$sr_comp_top10[8:nrow(shp)] <- FALSE

shp <- shp[order(shp$sr_added, decreasing=T), ]
shp$sr_comp_50 <- shp$sr_complementarity
shp$sr_comp_50[10:nrow(shp)] <- FALSE

# PD top 10 and comp 50%
shp <- shp[order(shp$pd_added, decreasing=T), ]
shp$pd_comp_50 <- shp$pd_complementarity
shp$pd_comp_50[24:nrow(shp)] <- FALSE

shp <- shp[order(shp$pd_added, decreasing=T), ]
shp$pd_comp_top10 <- shp$pd_complementarity
shp$pd_comp_top10[8:nrow(shp)] <- FALSE

# # PDE top 10 and comp 50%
shp <- shp[order(shp$PDE, decreasing=T), ]
shp$pde_comp_50 <- shp$pde_complementarity
shp$pde_comp_50[7:nrow(shp)] <- FALSE




# Top 2.5%  -----------------------------------------------------------
# complementarity vs top 2.5% values

# # Get total PD
# res = lapply(subphy, function(x){phyloregion::PD(submat, x)})
# # attach mean to shapefile
# tmp <- data.table(LEVEL3_COD = names(res[[1]]),
#                   PD_obs = rowMeans(as.data.table(res)))
# shp <- merge(shp, tmp, all.x=T)
# 
# # Get species richness
# tmp <- data.table(LEVEL3_COD = row.names(submat),
#                   richness = rowSums(submat))
# shp <- merge(shp, tmp, all.x=T)
# 
# shp$topsr = as.factor(phyloregion::hotspots(shp$richness))
# shp$toppd = as.factor(phyloregion::hotspots(shp$PD_obs))
# 
# View(st_drop_geometry(shp[,c("LEVEL3_COD", "LEVEL_NAME", "PD_obs", "richness",
#                              "topsr", "toppd", "sr_added", "pd_added")]))
# 


# save shp ------------------------------------------------------------
#saveRDS(shp, "data/pol_shp.rds")






# MAPS --------------------------------------------------------------------

shp = readRDS("data/pol_shp.rds")

# plot parameters
my_projection <- "+proj=wintri +datum=WGS84 +no_defs +over"

shp <- st_as_sf(shp)
shp <- st_transform(shp, crs=my_projection)
shp <- shp[!shp$LEVEL3_COD=="BOU",]

min.area <- 6e+9; min.area <- units::as_units(min.area, "m2")
shp$area <- st_area(shp)
thicc_lines <- shp[which(shp$area<min.area),]

bc <- c("#35abc4", "#4b9e31", "#eeea40")
my_projection <- "+proj=wintri +datum=WGS84 +no_defs +over"
# grat_wintri <-
#   sf::st_graticule(lat = c(-89.9, seq(-80, 80, 20), 89.9)) %>%
#   lwgeom::st_transform_proj(crs = my_projection)

grat_wintri <-
  st_graticule(lat = c(-98.9, 89.9),
               lon = c(-179.9, 179.9)) %>%
  st_transform_proj(crs = my_projection)

# set map theme
theme_set(theme_void()+
            theme(legend.position = c(0.1, 0.2),
                  legend.key.height = unit(6,"mm"),
                  legend.key.width = unit(4,"mm"),
                  legend.background = element_blank(),
                  legend.key = element_blank(),
                  legend.text.align=1,
                  panel.background = element_blank(),
                  panel.border = element_blank(),
                  text = element_text(size = 10)))

## FIG 2 ----------------------------------------------------------------------



(sr_top25 <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=topsr), lwd=0.25/.pt, col="gray95") + 
    coord_sf(expand=F, datum=NULL)+
    scale_fill_manual("top 2.5% SR", values=c(NA, "grey50"), na.value="grey80", labels=c('excluded', 'included'))+
    ggtitle("Political countries with top 2.5% SR") + 
    theme(legend.text.align=0, 
          legend.title=element_blank())
)

(pd_top25 <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=toppd), lwd=0.25/.pt, col="gray95") + 
    scale_fill_manual("top 2.5% PD", values=c(NA, bc[1]), na.value="grey80", labels=c('excluded', 'included'))+
    ggtitle("Political countries with top 2.5% PD")+
    coord_sf(expand=F, datum=NULL) + 
    theme(legend.text.align=0, 
          legend.title=element_blank())
)

(sr_comp10 <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=sr_comp_top10), lwd=0.25/.pt, col="gray95") + 
    scale_fill_manual("top 7 SR", values=c(NA, "grey50"), na.value="grey80", labels=c('excluded', 'included'))+
    coord_sf(expand=F, datum=NULL)+
    ggtitle("Top 7 most SR contributing political countries") + 
    theme(legend.text.align=0, 
          legend.title=element_blank())
)

(pd_comp10 <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=pd_comp_top10), lwd=0.25/.pt, col="gray95") + 
    scale_fill_manual("top 7 PD", values=c(NA, bc[1]), na.value="grey80", labels=c('excluded', 'included'))+
    coord_sf(expand=F, datum=NULL)+
    ggtitle("Top 7 most PD contributing political countries") + 
    theme(legend.text.align=0, 
          legend.title=element_blank())
)

plot_grid(
          sr_top25+theme(plot.title = element_text(hjust = 0.5)),
          pd_top25+theme(plot.title = element_text(hjust = 0.5)),
          sr_comp10+theme(plot.title = element_text(hjust = 0.5)),
          pd_comp10+theme(plot.title = element_text(hjust = 0.5)),
          ncol = 2, labels=c("A","B","C","D"), label_fontface=1, scale=1)
ggsave("figures/fig2_political_countries.png", width=10, height=6.5, units = "in", dpi = 300, bg = "white")


# plot_grid(sr_top25+theme(plot.title = element_text(hjust = 0.5)),
#           pd_top25+theme(plot.title = element_text(hjust = 0.5)),
#           sr_comp10+theme(plot.title = element_text(hjust = 0.5)),
#           pd_comp10+theme(plot.title = element_text(hjust = 0.5)),
#           ncol = 2, labels=c("A","B","C","D"), label_fontface=1, scale=1)
# ggsave("figures/fig2.png", width=10, height=6.5, units = "in", dpi = 300, bg = "white")


## FIG 3 ----------------------------------------------------------------------

(sr_comp_50_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=sr_comp_50), lwd=0.25/.pt, col="gray95") + 
    coord_sf(expand=F, datum=NULL)+
    scale_fill_manual("50% SR", values=c(NA, "grey50"), na.value="grey80", labels=c('excluded', 'included'))+
    ggtitle("50% global species richness") + 
   theme(legend.text.align=0, 
         legend.title=element_blank())
)

(pd_comp_50_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=pd_comp_50),lwd=0.25/.pt, col="gray95") + 
    coord_sf(expand=F, datum=NULL)+
    scale_fill_manual("50% PD", values=c(NA, bc[1]), na.value="grey80", labels=c('excluded', 'included'))+
    ggtitle("50% global PD") + 
    theme(legend.text.align=0, 
          legend.title=element_blank())
)

# (pde_comp_50_map <- ggplot(shp) + 
#     geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
#     geom_sf(aes(fill=pde_comp_50), lwd=0.25/.pt, col="gray95") + 
#     coord_sf(expand=F, datum=NULL)+
#     scale_fill_manual("50%PDE", values=c(NA, bc[2]), na.value="grey80")+
#     ggtitle("50% PD endemism with 12 bot countries")
# )

plot_grid(sr_comp_50_map+theme(plot.title = element_text(hjust = 0.5)),
          pd_comp_50_map+theme(plot.title = element_text(hjust = 0.5)),
          ncol=2, labels=c("A","B","C"), label_fontface=1, scale=1)
ggsave("figures/halflife_political.png", width=10, height=3.25, dpi=300, bg="white")








# STATS --------------------------------------------------------------------

# Global total SR or PD included in these countries

# SR top2.5% hotspots
sum(shp$sr_added_perc[shp$topsr==1])

# PD top2.5% hotspots
sum(shp$pd_added_perc[shp$toppd==1])


# Top 7 SR complementarity countries
sum(sort(shp$sr_added_perc, decreasing=T)[1:10])

# Top 7 PD complementarity countries
sum(sort(shp$pd_added_perc, decreasing=T)[1:10])



