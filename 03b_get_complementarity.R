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


# Load data

# only if still calculating things
#load("PD_nullmodel/comm_and_phy.RData")
source("99_functions.R")

shp <- readRDS("data/fin_shp.rds")




# Get SR complementarity -----------------------------------------------------

# greed <- sr_greedy(submat, m=0) #, n=10
# saveRDS(greed, "data/sr_complementarity.rds")

greed = readRDS("data/sr_complementarity.rds")
res <- data.table(LEVEL3_COD = names(greed[[1]]),
                  sr_complementarity = greed[[1]])

tmp = data.table(sr_added = greed[[2]],
                 LEVEL3_COD = greed[[3]])
res <- merge(res, tmp, all=T)
shp <- merge(shp, res, all.x=T)


# Get PD complementarity -----------------------------------------------------

# res <- list()
# for(i in 1:100){
#   greed <- pd_greedy(submat, phylo=subphy[[i]])
#   res[[i]] <- greed
#   cat(i, "\r")
# }
# saveRDS(res, "data/pd_complementarity.rds")
# gc()

# read results:
fnames <- dir("PDcomp")[grep("pd_complementarity", dir("PDcomp"))]
lpd <- list()
for(i in 1:100){
  lpd[[i]] <- readRDS(paste0("PDcomp/", fnames[i]))  
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

# merge into shapefile
shp <- merge(shp, res, all.x=T)



# Get PD endemism complementarity -------------------------------------------

# PDendemism is already its complementarity value, also the algorithm does not
# work since it alters absences of already covered countries, artificially
# increasing other countries endemism values. Use raw PDE values, ordered.

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
# To conserve 50% of SR, we need 15 countries

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
# To conserve 50% of PD, we need 33 countries

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
# To conserve 50% of PD, we need 12 countries


# scatterplot decrease added percent
bc <- c("#35abc4", "#4b9e31", "#eeea40")
plot(sort(shp$sr_added_perc, decreasing=T), type="b", ylim=c(0,0.1),main="black=SR, blue=PD, green=PDendemism")
points(sort(shp$pd_added_perc, decreasing=T), type="b", col=bc[1])
#points(sort(shp$pde_added_perc, decreasing=T), type="b", col=bc[2])

# order(shp$sr_added_perc, decreasing=T)
# shp2 = shp[order(shp$sr_added_perc, decreasing=T),]
# o1 = data.table(level = shp2$LEVEL3_COD, position = seq(nrow(shp2)))
# shp2 = shp[order(shp$pd_added_perc, decreasing=T),]
# o2 = data.table(level = shp2$LEVEL3_COD, position = seq(nrow(shp2)))
tmp = data.table(sr=sort(shp$sr_added_perc, decreasing=T), 
                 pd=na.omit(shp$pd_added_perc[order(shp$sr_added_perc, decreasing=T)]))


plot(tmp$sr, type="b", ylim=c(0,0.1),main="black=SR, blue=PD, green=PDendemism")
points(tmp$pd, type="b", col=bc[1])



# FIG 4 (barplot) --------------------------
pdt = data.table(percent_captured = seq(.1, 1, by=.1),
           SR = perc_sr, 
           PD = perc_pd, 
           PD_endemism = perc_pde)
pdt = melt(pdt, id="percent_captured")

(barplot = ggplot(pdt, aes(x= percent_captured, y=value, fill=variable, label=value))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_discrete("")+
  ylab("Botanical countries")+
  xlab("Proportion captured")+
  scale_fill_manual("", values=c("grey50", bc[1], bc[2]), labels=c('SR', 'PD', 'PD endemism'))+
  theme(legend.position=c(.25,.9), 
        legend.background=element_blank())
)
ggsave("figures/barplot_percentages_bot_countries.png", width=3.5, height=3.5, dpi=600, bg="white")
  






## Define top 10 / 50% -------------------------------------------------------

# SR top 10 and 50%
shp <- shp[order(shp$sr_added, decreasing=T), ]
shp$sr_comp_top10 <- shp$sr_complementarity
shp$sr_comp_top10[11:nrow(shp)] <- FALSE

shp <- shp[order(shp$sr_added, decreasing=T), ]
shp$sr_comp_50 <- shp$sr_complementarity
shp$sr_comp_50[16:nrow(shp)] <- FALSE

# PD top 10 and comp 50%
shp <- shp[order(shp$pd_added, decreasing=T), ]
shp$pd_comp_50 <- shp$pd_complementarity
shp$pd_comp_50[34:nrow(shp)] <- FALSE

shp <- shp[order(shp$pd_added, decreasing=T), ]
shp$pd_comp_top10 <- shp$pd_complementarity
shp$pd_comp_top10[11:nrow(shp)] <- FALSE

# PDE top 10 and comp 50%
shp <- shp[order(shp$PDE, decreasing=T), ]
shp$pde_comp_50 <- shp$pde_complementarity
shp$pde_comp_50[13:nrow(shp)] <- FALSE




# Top 2.5%  -----------------------------------------------------------
# complementarity vs top 2.5% values

shp$topsr = as.factor(phyloregion::hotspots(shp$richness))
shp$toppd = as.factor(phyloregion::hotspots(shp$PD_obs))


# _____________________ --------
# save shp ------------------------------------------------------------
saveRDS(shp, "data/fin_shp.rds")
# _____________________ --------




# MAPS --------------------------------------------------------------------

# plot parameters
my_projection <- "+proj=wintri +datum=WGS84 +no_defs +over"

shp <- st_as_sf(shp)
shp <- st_transform(shp, crs=my_projection)
shp <- shp[!shp$LEVEL3_COD=="BOU",]
#View(st_drop_geometry(shp[,c("LEVEL3_COD", "sr_complementarity")]))
min.area <- 6e+9; min.area <- units::as_units(min.area, "m2")
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
    ggtitle("Bot. countries with top 2.5% SR") + 
    theme(legend.text.align=0, 
          legend.title=element_blank())
)

(pd_top25 <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=toppd), lwd=0.25/.pt, col="gray95") + 
    scale_fill_manual("top 2.5% PD", values=c(NA, bc[1]), na.value="grey80", labels=c('excluded', 'included'))+
    ggtitle("Bot. countries with top 2.5% PD")+
    coord_sf(expand=F, datum=NULL) + 
    theme(legend.text.align=0, 
          legend.title=element_blank())
)

(sr_comp10 <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=sr_comp_top10), lwd=0.25/.pt, col="gray95") + 
    scale_fill_manual("top 10 SR", values=c(NA, "grey50"), na.value="grey80", labels=c('excluded', 'included'))+
    coord_sf(expand=F, datum=NULL)+
    ggtitle("Top 10 most SR contributing bot. countries")+ 
    theme(legend.text.align=0, 
          legend.title=element_blank())
)

(pd_comp10 <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=pd_comp_top10), lwd=0.25/.pt, col="gray95") + 
    scale_fill_manual("top 10 PD", values=c(NA, bc[1]), na.value="grey80", labels=c('excluded', 'included'))+
    ggtitle("Top 10 most PD contributing bot. countries")+
    coord_sf(expand=F, datum=NULL)+ 
    theme(legend.text.align=0, 
          legend.title=element_blank())
)

plot_grid(sr_top25+theme(plot.title = element_text(hjust = 0.5)),
          pd_top25+theme(plot.title = element_text(hjust = 0.5)),
          sr_comp10+theme(plot.title = element_text(hjust = 0.5)),
          pd_comp10+theme(plot.title = element_text(hjust = 0.5)),
          ncol = 2, labels=c("A","B","C","D"), label_fontface=1, scale=1)
ggsave("figures/fig2.png", width=10, height=6.5, units = "in", dpi = 300, bg = "white")



## FIG 3 ----------------------------------------------------------------------

(sr_comp_50_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=sr_comp_50), lwd=0.25/.pt, col="gray95") + 
    coord_sf(expand=F, datum=NULL)+
    scale_fill_manual("50%SR", values=c(NA, "grey50"), na.value="grey80", labels=c('excluded', 'included'))+
    ggtitle("50% SR in 15 bot. countries")+ 
   theme(legend.text.align=0, 
         legend.title=element_blank())
)

(pd_comp_50_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=pd_comp_50),lwd=0.25/.pt, col="gray95") + 
    coord_sf(expand=F, datum=NULL)+
    scale_fill_manual("50%PD", values=c(NA, bc[1]), na.value="grey80", labels=c('excluded', 'included'))+
    ggtitle("50% PD with 33 bot. countries")+ 
    theme(legend.text.align=0, 
          legend.title=element_blank())
)

(pde_comp_50_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=pde_comp_50), lwd=0.25/.pt, col="gray95") + 
    coord_sf(expand=F, datum=NULL)+
    scale_fill_manual("50%PDE", values=c(NA, bc[2]), na.value="grey80")+
    ggtitle("50% PD endemism with 12 bot. countries")
)

plot_grid(sr_comp_50_map+theme(plot.title = element_text(hjust = 0.5)),
          pd_comp_50_map+theme(plot.title = element_text(hjust = 0.5)),
          ncol=2, labels=c("A","B","C"), label_fontface=1, scale=1)
ggsave("figures/halflife.png", width=10, height=3.25, dpi=300, bg="white")







# _____________________ --------
# STATS --------------------------------------------------------------------

# Global total SR or PD included in these countries

# SR top2.5% hotspots
sum(shp$sr_added_perc[shp$topsr==1])

# PD top2.5% hotspots
sum(shp$pd_added_perc[shp$toppd==1])


# Top 10 SR complementarity countries
sum(sort(shp$sr_added_perc, decreasing=T)[1:10])

# Top 10 PD complementarity countries
sum(sort(shp$pd_added_perc, decreasing=T)[1:10])



# Total PD values for countries included in top 10 complementary: Would they be
# included in top 2.5%?




## TABLE 2 ----------
# table for hostpot traits. shows complementarity and top2.5% PD hotspots

tab = st_drop_geometry(shp[shp$pd_comp_top10 | shp$toppd==1,
                           c("LEVEL_NAME", "LEVEL3_COD", "topsr", "sr_comp_top10", 
                             "toppd","pd_comp_top10",
                             "richness", "PD_obs", "PDE")])
tab[,c("PD_obs", "PDE")] = round(tab[,c("PD_obs", "PDE")])
tab <- tab[order(tab$LEVEL_NAME),]
tab <- tab[order(tab$pd_comp_top10, decreasing=T),]
knitr::kable(tab, digits=2, format="simple", row.names=F)

# add numeric row for easy country ordering
tab$ord = c(1:nrow(tab))

# # top 2.5% PD table
# tab2 = st_drop_geometry(shp[shp$toppd==1,
#                            c("LEVEL_NAME", "LEVEL3_COD", "topsr", "sr_comp_top10", 
#                              "toppd","pd_comp_top10",
#                              "richness", "PD_obs", "PDE")])
# tab2[,c("PD_obs", "PDE")] = round(tab2[,c("PD_obs", "PDE")])
# tab2 <- tab2[order(tab2$LEVEL_NAME),]
# knitr::kable(tab2, digits=2, format="simple", row.names=F)

## Biomes ----
ggplot(tab)+
  geom_bar(aes(y=richness/max(shp$richness), x=factor(ord)), 
           stat="identity", fill="grey50", width=0.8)+
  scale_x_discrete(limits=rev)+
  coord_flip()+
  theme_classic()

# ggplot(tab2)+
#   geom_bar(aes(y=richness/max(shp$richness), x=LEVEL_NAME), 
#            stat="identity", fill="grey50", width=0.8)+
#   scale_x_discrete(limits=rev)+
#   coord_flip()+
#   theme_classic()

biomes <- readRDS("data/biomes_olson_ALL.rds")
biome_names <- c("(Sub)Tropical moist broadleaf forests", 
                 "(Sub)Tropical dry broadleaf forests", 
                 "(Sub)Tropical coniferous forests",
                 "Temperate broadleaf and mixed forests",
                 "Temperate coniferous forests", "Boreal forests/taiga",
                 "(Sub)Tropical grasslands, savannas, shrublands",
                 "Temperate grasslands, savannas, shrublands",
                 "Flooded_grasslands and savannas",
                 "Montane grasslands and shrublands", "Tundra",
                 "Mediterranean forests, woodlands, scrub",
                 "Deserts and xeric shrublands", "Mangroves", "Lakes", "Rock and ice")

names(biomes)[grepl("^X", names(biomes))] <- biome_names

# remove biomes that make up < 3% for clarity
for(i in 1:nrow(biomes)){
  tmp <- biomes[i,-1]
  tmp[which(tmp<0.03)] <- 0
  biomes[i,-1] <- tmp
}

rowSums(biomes[,-1], na.rm=T)
# scale to 100% to account for minor boundary inaccuracies
biomes[,-1] <- t(apply(biomes[,-1], 1, function(x){x/sum(x, na.rm=T)}) )
rowSums(biomes[,-1], na.rm=T)

hab.df <- tab[, c("LEVEL_NAME", "LEVEL3_COD", "ord")]
hab.df <- merge(hab.df, biomes, all.x=TRUE, by.x="LEVEL3_COD", by.y="country")
rowSums(hab.df[,-c(1,2,3)], na.rm=T) # should be all 1

hab.mlt <- melt(setDT(hab.df), id.vars=c("LEVEL_NAME", "LEVEL3_COD", "ord"))
hab.mlt <- hab.mlt[order(ord, value),]
hab.mlt$sort <- rep(length(unique(hab.mlt$variable)):1, length(unique(hab.mlt$LEVEL_NAME)))
# remove empty biomes
hab.mlt <- hab.mlt[hab.mlt$value!=0,]

hab.mlt <- droplevels(hab.mlt)
hab.mlt <- as.data.frame(hab.mlt)

# library(MetBrewer)
# mets = c("Archambault", "Cassatt1","Cassatt2" ,"Demuth","Derain","Egypt","Greek",
#          "Hiroshige","Hokusai2" ,"Hokusai3","Ingres","Isfahan1","Isfahan2","Java",
#          "Johnson","Kandinsky","Morgenstern","OKeeffe1","OKeeffe2","Pillement",
#          "Tam","Troy", "VanGogh3", "Veronese")

# my_pal <- c(scico::scico(9, alpha = .7, begin = 0.1, end = .9, direction = 1, palette = "batlow"), "#dbdbdb")
# my_pal0 = c("#0C335D", "#184E60", "#32665A", "#577646", "#808133", "#B28C32", "#E39856", "#FCA68C", "#FCB9C3", "#dbdbdb")
# my_pal2 <- c("#000000", "#FFFFFF", "#FFC0CB", "#1E90FF","#FFD700","#008000","#FFA500","#800080","#FF0000","#FFFF00")
my_pal = c("#1E3D14","#1F5B25","#669D62","#9CC184","#C2D6A4", "#B28C32", "#E39856", "#FCA68C", "#FCB9C3", "#0C335D")

bplot = ggplot(hab.mlt, aes(x=value, y=factor(ord), fill=variable)) +
  geom_bar(stat="identity", col="white", size=0.5)+
  scale_y_discrete(limits=rev)+
  scale_fill_manual("Biome", values=my_pal)+
#  scale_fill_manual(values=met.brewer(mets[24], 10))+
  guides(fill=guide_legend(nrow=4,byrow=TRUE))+
#  theme_classic()+
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"), 
        legend.spacing.y = unit(.1,"mm"),
        legend.text.align=0)
 #ggsave("figures/table1_biomes.pdf", width=6.5, height=5, units="in")

legend <- get_legend(bplot)                    
grid::grid.newpage()
grid::grid.draw(legend)

ggsave("figures/table1_biomes_legend.png", width=8, height=1, units="in")
#ggsave("figures/table1_biomes2.pdf", width=2, height=5, units="in")
bplot + theme(legend.position="none",
              plot.margin = unit(c(0, 0, 0, 0), "pt"))+
        coord_cartesian(expand=F)
ggsave("figures/table1_biomes2.png", width=2, height=3, units="in")


# hist for each biome type over all hotspots
tapply(hab.mlt$value, hab.mlt$variable, mean)
sort(round(tapply(hab.mlt$value, hab.mlt$variable, median),2), decreasing=T)
ggplot(hab.mlt, aes(x=value)) +
  geom_histogram()+
  facet_wrap("variable")+
  theme_classic()
  

  

# for top2.5%
hab.df <- tab2[, c("LEVEL_NAME", "LEVEL3_COD")]
hab.df <- merge(hab.df, biomes, all.x=TRUE, by.x="LEVEL3_COD", by.y="country")
rowSums(hab.df[,-c(1,2)], na.rm=T)

hab.mlt <- melt(setDT(hab.df))
hab.mlt <- hab.mlt[order(LEVEL_NAME, value),]
hab.mlt$sort <- rep(length(unique(hab.mlt$variable)):1, length(unique(hab.mlt$LEVEL_NAME)))
# remove empty biomes
hab.mlt <- hab.mlt[hab.mlt$value!=0,]

hab.mlt <- droplevels(hab.mlt)
hab.mlt <- as.data.frame(hab.mlt)

my_pal <- c(scico::scico(9, alpha = .7, begin = 0.1, end = .9, direction = 1, palette = "batlow"), "#dbdbdb")
ggplot(hab.mlt, aes(x=value, y=LEVEL_NAME, fill=variable)) +
  geom_bar(stat="identity", col="white", size=0.5)+
  scale_y_discrete(limits=rev)+
  scale_fill_manual("Biome", values=my_pal)+
  guides(fill=guide_legend(nrow=4,byrow=TRUE))+
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"), 
        legend.spacing.y = unit(.5,"mm"),
        legend.text.align=0)
#ggsave("figures/table1_biomes.pdf", width=6.5, height=5, units="in")
ggsave("figures/table1_biomes2.png", width=8, height=5, units="in")
#ggsave("figures/table1_biomes2.pdf", width=2, height=5, units="in")
ggsave("figures/table1_biomes22.png", width=2, height=5, units="in")


# hist for each biome type over all hotspots
sort(round(tapply(hab.mlt$value, hab.mlt$variable, median),2), decreasing=T)
ggplot(hab.mlt, aes(x=value)) +
  geom_histogram()+
  facet_wrap("variable")+
  theme_classic()




# ## BIOME MAP -----------------------------------------
# 
# # plot Olsons biomes
# 
# olson = vect("../TRFevol/data/olson_biomes/wwf_terr_ecos.shp")
# olson <- terra::project(olson, my_projection)
# olson <- st_as_sf(olson)
# 
# olson$BIOME[olson$BIOME==98] = 15
# olson$BIOME[olson$BIOME==98] = 16
# olson$BIOME = biome_names[olson$BIOME]
# ggplot() + 
#   geom_sf(data=olson, aes(fill=factor(BIOME)), lwd=0/.pt)+
#   coord_sf(expand=F, datum=NULL)
# 
# 
# pdhotspots = shp[shp$pd_comp_top10,]
# nopdhotspots = shp[!shp$pd_comp_top10,]
# 
# my_pal <- c(scico::scico(15, alpha = .7, begin = 0.1, end = .9, direction = 1, palette = "batlow"), "#dbdbdb")
# (biome_map <- ggplot() + 
#     geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
#     geom_sf(data=olson, aes(fill=factor(BIOME)), lwd=0/.pt)+
#     scale_fill_manual("Biome types", values=my_pal, na.value="grey80")+
#     ggtitle("PD complementary hotspots on Olson biomes")+
#     geom_sf(data=pdhotspots, fill=NA, col="black", lwd=1.5/.pt)+
#     coord_sf(expand=F, datum=NULL)+
#     guides(fill=guide_legend(ncol=1,byrow=F))+
#     theme(legend.text.align=0)
# )
# ggsave("figures/biome_map.png", width=8, height=5.25, units = "in", dpi = 300, bg = "white")
# 
# 
# 
# # cropped version
# b2 = biomes[biomes$country %in% pdhotspots$LEVEL3_COD,]
# b2 = b2!=0
# # delete those
# tmp = apply(b2[,-1], 1, function(x){which(x)})
# tmp = unlist(tmp)
# included_biomes = unique(gsub("[0-9]{2,3}\\.", "", names(tmp)))
# olson <- olson[olson$BIOME %in% included_biomes,]
# 
# my_pal <- c(scico::scico(12, alpha = .9, begin = 0.1, end = .9, direction = 1, palette = "batlow"))
# my_pal[7] = NA
# (biome_map <- ggplot() + 
#     geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
#     geom_sf(data=olson, aes(fill=factor(BIOME)), lwd=0/.pt)+
#     scale_fill_discrete("Biome types")+#, values=my_pal, na.value="grey80")+
#     ggtitle("PD complementary hotspots with Olson biomes")+
#     geom_sf(data=nopdhotspots, fill="grey80", lwd=0/.pt)+
#     geom_sf(data=pdhotspots, fill=NA, col="black", lwd=1.5/.pt)+
#     coord_sf(datum=NULL, xlim=c(-9390315, 14000000), ylim=c(-4500000, 4767000))+
#     guides(fill=guide_legend(ncol=3,byrow=T))+
#     theme(legend.text.align=0,
#           legend.position="bottom")
# )
# 
# (biome_map <- ggplot() + 
#     geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
#     geom_sf(data=olson, aes(fill=factor(BIOME)), lwd=0/.pt)+
#     scale_fill_manual("Biome types", values=my_pal, na.value="grey80")+
#     ggtitle("PD complementary hotspots on Olson biomes")+
#     geom_sf(data=pdhotspots, fill=NA, col="black", lwd=1.5/.pt)+
#     coord_sf(datum=NULL, xlim=c(-15090315, 15000000), ylim=c(-6500000, 5567000))+
#     guides(fill=guide_legend(ncol=1,byrow=F))+
#     theme(legend.text.align=0)
# )
# ggsave("figures/biome_map2.png", width=8, height=5.25, units = "in", dpi = 300, bg = "white")
# 
# #ggsave("figures/biome_map_cropped.png", width=8, height=5.25, units = "in", dpi = 300, bg = "white")

# # in panels???
# countries = pdhotspots$LEVEL3_COD
# nopdhotspot_names = paste0("nohotspots", countries)
# for(i in 1:length(countries)){
#   tmp = st_crop(olson, pdhotspots[i,])
#   assign(countries[i], tmp)
#   assign(nopdhotspot_names[i], st_crop(nopdhotspots, tmp))
# }
# 
# maps = paste0(countries, "_map") 
# for(index in 1:length(countries)){
# tmp = ggplot() + 
#     geom_sf(data=get(countries[index]), aes(fill=factor(BIOME)), lwd=0/.pt)+
#     scale_fill_manual("Biome types", values=my_pal, na.value="grey80")+
#     ggtitle(countries[index])+
#     geom_sf(data=get(nopdhotspot_names[index]), fill="grey80", lwd=0/.pt)+
#     geom_sf(data=pdhotspots[index,], fill=NA, col="black", lwd=1.5/.pt)+
#     coord_sf(expand=F, datum=NULL)+
#     guides(fill=guide_legend(ncol=1,byrow=F))+
#     theme(legend.text.align=0)
# assign(maps[index], tmp)
# }
# 
# mget(maps)
# plot_grid(get(maps[1]), get(maps[2]), get(maps[3]), get(maps[4]), get(maps[5]),
#           get(maps[6]), get(maps[7]), get(maps[8]), get(maps[9]), get(maps[10]),
#           nrow=2)




## Biome stats ----------------------------------------------------------------
tmp = st_drop_geometry(shp[,c("LEVEL3_COD", "pd_comp_top10")])
setDT(tmp)
tmp <- merge(tmp, biomes, all.x=TRUE, by.x="LEVEL3_COD", by.y="country")
tmp2 = melt(tmp, id.vars=c("LEVEL3_COD", "pd_comp_top10"), value.name = 'coverage', variable.name='biome')

(comp_biome = ggplot(tmp2)+
  geom_boxplot(aes(y=biome, x=coverage, fill=pd_comp_top10), outlier.alpha = .2)+
  scale_fill_manual("top PD\ncomp", values=c("white", "grey80"))+
  theme_classic()+
  ylab("")+
  theme(legend.position=c(.8,.9))
  #scale_y_discrete(position="right")
)
# single tests
# for(i in 1:length(unique(tmp2$biome))){
#   s = tmp2[tmp2$biome==unique(tmp2$biome)[i],]
#   if(kruskal.test(s$coverage, s$toppd)$p.value<0.05){
#     print(s$biome[1], max.levels=0)
#     print(kruskal.test(s$coverage, s$pd_comp_top10))
#     print(tapply(s$coverage, s$pd_comp_top10, psych::describe))
#   }
# }

# Top 2.5%
tmp = st_drop_geometry(shp[,c("LEVEL3_COD", "toppd")])
setDT(tmp)
tmp <- merge(tmp, biomes, all.x=TRUE, by.x="LEVEL3_COD", by.y="country")
tmp2 = melt(tmp, id.vars=c("LEVEL3_COD", "toppd"), value.name = 'coverage', variable.name='biome')

(top_biome = ggplot(tmp2)+
  geom_boxplot(aes(y=biome, x=coverage, fill=toppd), outlier.alpha = .2)+
  scale_fill_manual("top 2.5%\nPD", values=c("white", "grey80"))+
  theme_classic()+
  ylab("")+
  theme(legend.position=c(.8,.9),
    axis.text.y=element_blank())
)

plot_grid(comp_biome, top_biome, ncol=2, rel_widths=c(.63,.37),
          labels=c("A","B"), label_fontface=1, scale=1)

ggsave('figures/biome_proportions.png', width=9, height=5, dpi=300)

# single tests
for(i in 1:length(unique(tmp2$biome))){
  s = tmp2[tmp2$biome==unique(tmp2$biome)[i],]
  if(kruskal.test(s$coverage, s$toppd)$p.value<0.05){
    print(s$biome[1], max.levels=0)
    print(kruskal.test(s$coverage, s$toppd))
    print(tapply(s$coverage, s$toppd, psych::describe))
  }
}



### number of biomes correlated with total PD--------------------------------

# number of biomes per country
biomes <- readRDS("data/biomes_olson_ALL.rds")
biome_names <- c("(Sub)Tropical moist broadleaf forests", 
                 "(Sub)Tropical dry broadleaf forests", 
                 "(Sub)Tropical coniferous forests",
                 "Temperate broadleaf and mixed forests",
                 "Temperate coniferous forests", "Boreal forests/taiga",
                 "(Sub)Tropical grasslands, savannas, shrublands",
                 "Temperate grasslands, savannas, shrublands",
                 "Flooded_grasslands and savannas",
                 "Montane grasslands and shrublands", "Tundra",
                 "Mediterranean forests, woodlands, scrub",
                 "Deserts and xeric shrublands", "Mangroves", "Lakes", "Rock and ice")

names(biomes)[grepl("^X", names(biomes))] <- biome_names
b2 = biomes[,-1]
b2 = b2!=0
tmp = apply(b2[,-1], 1, function(x){length(which(x))})
tmpdf = data.table(LEVEL3_COD = biomes[,1],
                   number_biomes = tmp)
shp = merge(shp, tmpdf, all.x=T)

ggplot(shp, aes(x=number_biomes, PD_obs))+
  geom_point()+
  ylab("PD")+xlab("number of biomes")+
  geom_smooth()+
  theme_classic()
#ggsave("figures/number_biomes_PD.png", width=3.5, height=3.5, dpi=600, bg="white")

cor.test(shp$number_biomes, shp$PD_obs, method="s") # rho = 0.40



# complementarity hotspots have more biome types than: 
  # a) non-complementarity hotspots
  tmp <- shp[shp$pd_comp_top10 | shp$toppd==1,]
  # ggplot(tmp)+
  #   geom_boxplot(aes(x=pd_comp_top10, y=mrd, fill=pd_comp_top10), show.legend = F,
  #                varwidth=T, width=1)+
  #   labs(y="Net diversification rate (mean root distance)",
  #        x="PD top 10")+
  #   scale_fill_manual(values=c("grey", bc))+
  #   theme_classic()
  kruskal.test(tmp$number_biomes, tmp$pd_comp_top10)
  kruskal.test(tmp$number_biomes, tmp$toppd)
  tapply(tmp$number_biomes, tmp$pd_comp_top10, psych::describe)
  tapply(tmp$number_biomes, tmp$toppd, psych::describe)

  
  # b) all other countries
  kruskal.test(shp$number_biomes, shp$pd_comp_top10)
  tapply(shp$number_biomes, shp$pd_comp_top10, psych::describe)

# PD total hotspots have more biome types than...
  # b) all other countries
  
  kruskal.test(shp$number_biomes, shp$toppd)
  tapply(shp$number_biomes, shp$toppd, psych::describe)




## Div rates -------------------------------------------------------------------

kruskal.test(shp$mrd, shp$pd_comp_top10)
kruskal.test(shp$mrd, shp$pd_comp_50)
kruskal.test(shp$mrd, shp$toppd)

#kruskal.test(shp$mrd, shp$sr_comp_top10)
#kruskal.test(shp$mrd, shp$sr_comp_50)


bc <- c("#35abc4", "#4b9e31", "#eeea40")
bp1 = ggplot(shp)+
  geom_boxplot(aes(x=pd_comp_top10, y=mrd, fill=pd_comp_top10), show.legend = F,
               varwidth=T, width=1)+
  labs(y="",
       x="PD top 10")+
  scale_fill_manual(values=c("grey", bc))+
  #annotate("text", x=1.5, y=120.7, label="P=0.038", size=2)+
  #scale_x_discrete(labels=c("50%", "top 10 PD"))+
  #geom_segment(aes(x=1,y=120, xend=2,yend=120), lwd=0.1) #, arrow = arrow(length = unit(0.1,"cm"))+
  theme_classic()

bp2 = ggplot(shp)+
  geom_boxplot(aes(x=pd_comp_50, y=mrd, fill=pd_comp_50), show.legend = F,
               varwidth=T, width=1)+
  labs(y="",
       x="PD 50%")+
  scale_fill_manual(values=c("grey", bc))+
  #annotate("text", x=1.5, y=120.9, label="P=0.027", size=4)+
  #scale_x_discrete(labels=c("50%", "top 10 PD"))+
  #geom_segment(aes(x=1,y=120, xend=2,yend=120), lwd=0.1)+ #, arrow = arrow(length = unit(0.1,"cm"))+
  theme_classic()

bp3 = ggplot(shp)+
  geom_boxplot(aes(x=toppd, y=mrd, fill=toppd), show.legend = F,
               varwidth=T, width=1)+
  labs(y="Net diversification rate (mean root distance)",
       x="top 2.5% PD")+
  #annotate("text", x=1.5, y=120.9, label="P=0.04", size=4)+
  #scale_x_discrete(labels=c("50%", "top 10 PD"))+
  #geom_segment(aes(x=1,y=120, xend=2,yend=120), lwd=0.1)+ #, arrow = arrow(length = unit(0.1,"cm"))+
  theme_classic()


plot_grid(bp3, bp1, bp2, ncol=3)
ggsave("figures/hotspots_boxplots_MRD.png", width=8, height=4, units = "in", dpi = 300, bg="white")



## Hotspots VS non-hotspots - threats ---------------------------

shp2.tmp <- st_drop_geometry(shp[,c("pd_comp_top10", "deforestation2", "hfp.1", "mat_change", "pre_change")])
shp2.tmp[,-1] <- apply(shp2.tmp[,-1], 2, normalized)
shp2.tmp <- data.table::melt(shp2.tmp)

bat <- scico::scico(palette="hawaii", n=4, begin=0, alpha=1)
bat2 <- c(bat[3], bat[2], bat[1], bat[4])
ggplot(shp2.tmp, aes(x=pd_comp_top10, y=value, fill=variable))+
  geom_boxplot(varwidth=F)+
  scale_fill_manual("Threat type", values=bat2, labels=c("deforestation", "HFP", "MAT change", "PRE change"))+
  scale_y_continuous("Scaled threat values", lim=c(0,1.03))+
  xlab("")+
  theme_classic()+
  theme(legend.position=c(.8,.89), 
        legend.background=element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE, keyheight=unit(6,"mm")))+
  annotate("text", x=1.2,y=1.03, label="italic(P) == 0.036", size=4, parse=T)+
  geom_segment(aes(x=0.7,y=1.014, xend=1.7,yend=1.015), lwd=0.3)+
  scale_x_discrete(labels=c("no hotspot", "PD comp top 10"))
ggsave("figures/threat_boxplots.png", width=7, height=6, units = "in", dpi = 300, bg="white")

kruskal.test(shp$hfp.1, shp$pd_comp_top10) #nope
kruskal.test(shp$mat_change, shp$pd_comp_top10) # nope
kruskal.test(shp$pre_change, shp$pd_comp_top10) # nope
kruskal.test(shp$deforestation2, shp$pd_comp_top10) # yes

# quantify
tapply(shp$deforestation2, shp$pd_comp_top10, psych::describe)






## AREA bias? ----------------

#View(st_drop_geometry(shp[,c("LEVEL_NAME", "area", "pd_comp_top10", "pd_comp_50")]))
#shp$LEVEL_NAME[shp$sr_comp_top10]
shp$LEVEL_NAME[shp$pd_comp_top10]
#shp$LEVEL_NAME[shp$sr_comp_50]
shp$LEVEL_NAME[shp$pd_comp_50]

# comp top 10
kruskal.test(as.numeric(shp$area), shp$pd_comp_top10)
p1 = ggplot(shp, aes(x=pd_comp_top10, y=as.numeric(area), fill=pd_comp_top10))+
  geom_boxplot(show.legend=F)+
  scale_fill_manual(values=c("grey", bc[1]))+
  scale_y_log10()+
  theme_classic()+
  ylab("Area")

kruskal.test(shp$area, shp$pd_comp_50)
p2= ggplot(shp, aes(x=pd_comp_50, y=as.numeric(area), fill=pd_comp_50))+
  geom_boxplot(show.legend=F)+
  scale_fill_manual(values=c("grey", bc[1]))+
  scale_y_log10()+
  theme_classic()+
  ylab("Area")

plot_grid(p1,p2, ncol=2)
ggsave("figures/area_boxplots.png", width=6, height=4, units = "in", dpi = 300, bg="white")








# _____________________ --------
# Conservation hotspot maps --------------------------------------------------

ha_united <- st_read("data/hotspots_merged_cleaned.gpkg")
ha_united <- st_wrap_dateline(ha_united, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
ha_united <- st_transform(ha_united, my_projection)
ha_united3 <- st_simplify(ha_united)

(pd_comp_hotspot_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=pd_comp_top10),lwd=0.25/.pt, col="gray95") + 
#    geom_sf(data=thicc_lines, lwd=1, aes(col=pd_comp_top10), show.legend=F)+
    theme_void()+
    scale_fill_manual(values=c(NA, bc[1]), na.value="grey80")+
#    scale_color_manual(values=c(NA, bc[1]), na.value="grey80")+
    theme(legend.position = "c(0.2, 0.3)",
          legend.key.height = unit(6,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text.align=1,
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )+
  geom_sf(data=ha_united3, col=NA, fill=bc[3], alpha=.3, size=.5)+ #alpha("white",0.9)
  coord_sf(expand=F, datum=NULL)
)
#ggsave("figures/PDhalflife_HS_map.png", width=6, height=3.5, units = "in", dpi = 600, bg = "white")

(sr_comp_hotspot_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=sr_comp_top10),lwd=0.25/.pt, col="gray95") + 
    #    geom_sf(data=thicc_lines, lwd=1, aes(col=sr_comp_top10), show.legend=F)+
    theme_void()+
    scale_fill_manual(values=c(NA, bc[1]), na.value="grey80")+
    #    scale_color_manual(values=c(NA, bc[1]), na.value="grey80")+
    theme(legend.position = "c(0.2, 0.3)",
          legend.key.height = unit(6,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text.align=1,
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )+
    geom_sf(data=ha_united3, col=NA, fill=bc[3], alpha=.3, size=.5)+ #alpha("white",0.9)
    coord_sf(expand=F, datum=NULL)
)

(pde_comp_hotspot_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=pde_complementarity2),lwd=0.25/.pt, col="gray95") + 
    #    geom_sf(data=thicc_lines, lwd=1, aes(col=pde_complementarity2), show.legend=F)+
    theme_void()+
    scale_fill_manual(values=c(NA, bc[1]), na.value="grey80")+
    #    scale_color_manual(values=c(NA, bc[1]), na.value="grey80")+
    theme(legend.position = "c(0.2, 0.3)",
          legend.key.height = unit(6,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text.align=1,
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )+
    geom_sf(data=ha_united3, col=NA, fill=bc[3], alpha=.3, size=.5)+ #alpha("white",0.9)
    coord_sf(expand=F, datum=NULL)
)




# PD top2.5% vs conservation hotspots 
(pd_top25_hotspot_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=toppd),lwd=0.25/.pt, col="gray95") + 
    #    geom_sf(data=thicc_lines, lwd=1, aes(col=pd_comp_top10), show.legend=F)+
    theme_void()+
    scale_fill_manual(values=c(NA, bc[1]), na.value="grey80")+
    #    scale_color_manual(values=c(NA, bc[1]), na.value="grey80")+
    theme(legend.position = "c(0.2, 0.3)",
          legend.key.height = unit(6,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text.align=1,
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )+
    geom_sf(data=ha_united3, col=NA, fill=bc[3], alpha=.3, size=.5)+ #alpha("white",0.9)
    coord_sf(expand=F, datum=NULL)
)
#ggsave("figures/PDtop25_HS_map.png", width=6, height=3.5, units = "in", dpi = 600, bg = "white")


(sr_top25_hotspot_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray60", size = 0.25/.pt) + 
    geom_sf(aes(fill=toppd),lwd=0.25/.pt, col="gray95") + 
    #    geom_sf(data=thicc_lines, lwd=1, aes(col=pd_comp_top10), show.legend=F)+
    theme_void()+
    scale_fill_manual(values=c(NA, bc[1]), na.value="grey80")+
    #    scale_color_manual(values=c(NA, bc[1]), na.value="grey80")+
    theme(legend.position = "c(0.2, 0.3)",
          legend.key.height = unit(6,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text.align=1,
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )+
    geom_sf(data=ha_united3, col=NA, fill=bc[3], alpha=.3, size=.5)+ #alpha("white",0.9)
    coord_sf(expand=F, datum=NULL)
)












# old stuff -------------------------------------------------------------------


# alternative: Get average complementarity compared to a random set 
# 
# added_main <- c()
# added_sd <- c()
# set.seed(1354)
# setsize <- 50
# for(i in 1:nrow(submat)){
#   cell_sr <- names(which(submat[i,]>0))
#   added <- c()
#   # random sample 100 times
#   for(j in 1:100){
#     # get random SR for set of x countries
#     mat <- submat[sample(c(1:nrow(submat)), setsize), ]
#     # count all columns with at least one entry
#     species_random <- names(which(colSums(mat)>0))
#     # how many cell species are unique to the random sample?
#     un <- cell_sr[which(!cell_sr %in%  species_random)]
#     added <- c(added, length(un))
#     #cat(j, "\r")
#   }
#   added_main <- c(added_main, mean(added))
#   added_sd <- c(added_sd, sd(added))
#   
#   cat(i, "\r")
# }
# 
# #saveRDS(list(added_main, added_sd), "data/added_sr.rds")
# 
# tmp <- data.table(LEVEL3_COD = row.names(submat), 
#                   added_main, 
#                   added_sd)
# 
# shp <- merge(shp, tmp, all.x=T)
# ggplot(shp)+
#   geom_sf(aes(fill=added_main),lwd=0, col=NA) + 
#   #  geom_sf(data=thicc_lines, lwd=1, aes(col=pd_added), show.legend=F)+
#   coord_sf(expand=F, datum=NULL)+
#   theme_void()+
#   theme(legend.position = c(0.2, 0.3))+
#   ggtitle("Average number species added to a random sample of 50 countries")
# plot(shp$added_main, shp$added_sd)
# # Get WE complementarity -----------------------------------------------------
# 
# we_greedy <- function(species_matrix, m=0, n = nrow(species_matrix)) {
#   
#   # set to store the cells
#   cell_set <- logical(nrow(species_matrix))
#   names(cell_set) <- row.names(species_matrix)
#   # store the number of species in each cell
#   species_count <- phyloregion::weighted_endemism(species_matrix)
#   # store species count in each iteration
#   species <- c()
#   area <- c()
#   
#   # Loop until all species are represented (species count=0)
#   while (sum(species_count) > m & sum(cell_set) < n) {
#     # Find the cell with the maximum number of species and add it to set
#     best_cell <- which.max(species_count)
#     cell_set[best_cell] <- TRUE
#     # store species count with country
#     species <- c(species, max(species_count))
#     area <- c(area, names(best_cell))
#     # Update the species count for the remaining cells
#     ## set species (=columns) that are represented in best cell(=row) to 0
#     counted_species <- which(species_matrix[best_cell, ]==1)
#     species_matrix[, counted_species] <- 0
#     
#     species_count <-  phyloregion::weighted_endemism(species_matrix)
#   }
#   return(list(cell_set, species, area))
# }
# 
# greed <- sr_greedy(submat, m=0) #, n=10 
# Get PD complementarity corrected for AREA ----------------------------------------
area_dt = data.table(bc = shp$LEVEL3_COD, area=shp$area)
# 
# # testing area
# species_matrix = submat
# phylo=subphy[[1]]
# area_dt = area_dt

pd_area_greedy <- function(species_matrix, phylo, m=0, n = nrow(species_matrix),
                           area_dt) { # tmp = area and botcountry dt
  
  # set to store the cells
  cell_set <- logical(nrow(species_matrix))
  names(cell_set) <- row.names(species_matrix)
  # store the PD in each cell
  pd <- phyloregion::PD(species_matrix, phylo)
  # store species count in each iteration
  pd_number <- c()
  botcountry <- c()
  
  # Loop until all PD is represented (sum(pd)=0)
  while (sum(pd) > m & sum(cell_set) < n) {
    
    # Find the cell with the maximum total PD PER AREA and add it to set
    tmp = merge(area_dt, data.table(bc=names(pd), pd), all.y=T) # MERGE, keep all to not mess up order
    pd_area = tmp$pd/tmp$area # PD per country per area
    
    # pick the best relative PD
    best_cell <- which.max(pd_area)[1]
    cell_set[best_cell] <- TRUE
    
    # store pd with country
    pd_number <- c(pd_number, pd[best_cell]) # store absolute PD, not relative PD
    botcountry <- c(botcountry, names(cell_set[best_cell]))
    
    # Update the pd for the remaining cells
    ## set species (=columns) that are represented in best cell(=row) to absent. 
    counted_species <- which(species_matrix[best_cell, ]==1) # too many replacement errors
    ## add switch for empty best cells....
    if(length(counted_species)>0){
      species_matrix[, counted_species] <- 0
      pd <-  phyloregion::PD(species_matrix, phylo)
      message(paste(names(cell_set[best_cell]), " : ", max(pd))) # returns max pd to show progress
    }else{
      break
    }
  }
  return(list(cell_set, pd_number, botcountry))
}


#greed <- pd_area_greedy(submat, phylo=subphy[[1]], area_dt=area_dt)
#saveRDS(greed, "data/PD_area_comp.rds")
greed = readRDS("data/PD_area_comp.rds")

# results
res <- data.table(LEVEL3_COD = names(greed[[1]]),
                  pdarea_complementarity = greed[[1]])
tmp = data.table(pdarea_added = greed[[2]],
                 LEVEL3_COD = greed[[3]],
                 add_order = seq(1:length(greed[[3]])))
res <- merge(res, tmp, all=T)

shp <- merge(shp, res, all.x=T)
thicc_lines <- shp[which(shp$area<min.area),]

View(st_drop_geometry(shp[,c("LEVEL_NAME", "pd_added", "pdarea_added")]))

# get percentages
shp$pda_added_perc = shp$pdarea_added/sum(shp$pdarea_added, na.rm=T)
perc_pda = c()
for(i in seq(.1, 1, by=.1)){
  shppda = na.omit(shp$pda_added_perc)
  tmp = c()
  while(sum(tmp) <= i & length(shppda)>0) {
    # take percentage and add to set
    tmp <- c(tmp, max(shppda, na.rm=T))
    shppda = shppda[-which.max(shppda)]
  }
  perc_pda = c(perc_pda, length(tmp))
  cat(i, "\r")
}
perc_pda

# plot
# cannot simply order for area added since this differs, order for add order
shp <- shp[order(shp$add_order, decreasing=F), ]
shp$pdarea_complementarity2 <- shp$pdarea_complementarity
shp$pdarea_complementarity2[57:nrow(shp)] <- FALSE

thicc_lines <- shp[which(shp$area<min.area),]
(pdarea_complementarity_map <- ggplot(shp) + 
    geom_sf(aes(fill=pdarea_complementarity2),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=pdarea_complementarity2), show.legend=F)+
    theme_void()+
    coord_sf(expand=F, datum=NULL)+
    theme(legend.position = c(0.2, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text.align=1,
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )+
    ggtitle("50% PD, contributors selected greedy but for max PD/area")
)
plot(shp$pdarea_added)

# make some cumulative plot
#plot(ecdf(shp$sr_added))
# ggplot(shp, aes(x = sr_added)) +
#   stat_function(fun = pnorm)
# 
# plot(abs(diff(sort(shp$sr_added_perc, decreasing=T))),
#      ylab="difference in number of added species")




# # Get SR complementarity with area guidance -----------------------------
area_dt = data.table(bc = shp$LEVEL3_COD, area=shp$area)
# 
# # testing area
# species_matrix = submat
# phylo=subphy[[1]]
# area_dt = area_dt

sr_area_greedy <- function(species_matrix, phylo, m=0, n = nrow(species_matrix),
                           area_dt) { # tmp = area and botcountry dt
  
  # set to store the cells
  cell_set <- logical(nrow(species_matrix))
  names(cell_set) <- row.names(species_matrix)
  # store the PD in each cell
  pd <- rowSums(species_matrix, phylo)
  # store species count in each iteration
  pd_number <- c()
  botcountry <- c()
  
  # STOP HERE
  
  # Loop until all PD is represented (sum(pd)=0)
  while (sum(pd) > m & sum(cell_set) < n) {
    
    # Find the cell with the maximum total PD PER AREA and add it to set
    tmp = merge(area_dt, data.table(bc=names(pd), pd), all.y=T) # MERGE, keep all to not mess up order
    pd_area = tmp$pd/tmp$area # PD per country per area
    
    # pick the best relative PD
    best_cell <- which.max(pd_area)[1]
    cell_set[best_cell] <- TRUE
    
    # store pd with country
    pd_number <- c(pd_number, pd[best_cell]) # store absolute PD, not relative PD
    botcountry <- c(botcountry, names(cell_set[best_cell]))
    
    # Update the pd for the remaining cells
    ## set species (=columns) that are represented in best cell(=row) to absent. 
    counted_species <- which(species_matrix[best_cell, ]==1) # too many replacement errors
    ## add switch for empty best cells....
    if(length(counted_species)>0){
      species_matrix[, counted_species] <- 0
      pd <-  phyloregion::PD(species_matrix, phylo)
      message(paste(names(cell_set[best_cell]), " : ", max(pd))) # returns max pd to show progress
    }else{
      break
    }
  }
  return(list(cell_set, pd_number, botcountry))
}


greed <- pd_area_greedy(submat, phylo=subphy[[1]], area_dt=area_dt)



# OpenChatAI

# constrained optimization

# load required packages
library(gapminder)
library(dplyr)

# group data by country and summarize area and temperature
country_summaries <- st_drop_geometry(shp[,c("LEVEL3_COD","richness", "area")])
  # gapminder %>% 
  # group_by(country) %>% 
  # summarize(total_area = sum(area),
  #           total_temp = sum(temp))

# define function to minimize area and maximize richness
obj_fun <- function(x) {
  area <- sum(country_summaries$area[x])
  richness <- sum(country_summaries$richness[x])
  return(list(area = area, richness = richness))
}

# set optimization constraints
n_countries <- 1
constraint_fun <- function(x) length(x) == n_countries

# perform optimization
result <- optimize(obj_fun, 
                   lower = rep(TRUE, nrow(country_summaries)), 
                   upper = rep(FALSE, nrow(country_summaries)), 
                   constraints = list("eq" = constraint_fun), 
                   maximum = TRUE)

# extract the countries that minimize area and maximize temperature
country_names <- country_summaries$country[result$ix]

# load required packages
library(gapminder)
library(dplyr)

# group data by country and summarize area and temperature
country_summaries <- gapminder %>% 
  group_by(country) %>% 
  summarize(total_pop = sum(pop),
            total_lifeExp = sum(lifeExp))

# define function to minimize area and maximize temperature
obj_fun <- function(x) {
  pop <- sum(country_summaries$total_pop[x])
  return(list(pop = pop, life = life))
}

# set optimization constraints
n_countries <- 5
constraint_fun <- function(x) length(x) == n_countries

# perform optimization
result <- optimize(obj_fun, 
                   lower = 4, 
                   upper = 148, 
                   constraints = list("eq" = constraint_fun), 
                   maximum = TRUE)

# extract the countries that minimize area and maximize temperature
country_names <- country_summaries$country[result$ix]



# load GA package
library(GA)

# define objective function
obj_fun <- function(x, y) { # takes two variables
  # extract variables
  x1 <- x
  x2 <- y
  # define weights
  w1 <- 0.5 # weight for x
  w2 <- 0.5 # weight for y
  # calculate weighted sum
  obj_val <- w1 * x1 - w2 * x2
  # return objective value
  return(obj_val)
}

# define parameter bounds and type
param_bounds <- matrix(c(0, 10), ncol = 2, byrow = TRUE)
param_type <- rep("real", 2)

# perform optimization
result <- ga(type = "real-valued", 
             fitness = obj_fun, 
             lower = -10, 
             upper = 10,
             x=shp$area, y=shp$richness)

# extract optimal solution
opt_sol <- result$population[1, ]



f <- function(x)  abs(x)+cos(x)
curve(f, -20, 20)

fitness <- function(x) -f(x)
GA <- ga(type = "real-valued", fitness = fitness, lower = -100, upper = 100)
summary(GA)
plot(GA)



shp2 = shp[!shp$LEVEL3_COD=="ANT",]
par(mfrow=c(1,2))

plot(sqrt(as.numeric(shp2$area)), shp2$PD_obs, pch=20, 
     col=c("black","red")[as.numeric(shp2$pd_comp_top10)+1])

plot(sqrt(as.numeric(shp2$area)), shp2$PD_obs, pch=20,
    col=c("black","red")[as.numeric(shp2$pd_comp_50)+1])


plot(sqrt(as.numeric(shp2$area)), shp2$pd_added, pch=20, 
     col=c("black","red")[as.numeric(shp2$pd_comp_top10)+1])
plot(sqrt(as.numeric(shp2$area)), shp2$pd_added, pch=20, 
     col=c("black","red")[as.numeric(shp2$pd_comp_50)+1])



# (sr_complementarity_bi <- ggplot(shp) + 
#     geom_sf(aes(fill=sr_added),lwd=0, col=NA) + 
#     geom_sf(data=thicc_lines, lwd=1, aes(col=sr_added), show.legend=F)+
#     theme_void()+
#     scale_fill_continuous(trans="sqrt")+
#     scale_color_continuous(trans="sqrt")+
#     coord_sf(expand=F, datum=NULL)+
#     theme(legend.position = c(0.2, 0.3),
#           legend.key.height = unit(6,"mm"),
#           legend.key.width = unit(4,"mm"),
#           legend.background = element_blank(),
#           legend.key = element_blank(),
#           legend.text.align=1,
#           panel.background = element_blank(),
#           panel.border = element_blank(),
#           text = element_text(size = 10)
#     )+
#     ggtitle("number of new species added to the set, starting with CLM")
# )

# (pd_complementarity_bi <- ggplot(shp) + 
#     geom_sf(aes(fill=pd_added),lwd=0, col=NA) + 
#     geom_sf(data=thicc_lines, lwd=1, aes(col=pd_added), show.legend=F)+
#     scale_fill_continuous(trans="sqrt")+
#     scale_color_continuous(trans="sqrt")+
#     coord_sf(expand=F, datum=NULL)+
#     ggtitle("amount of PD added to the set, starting with CLM, average for 100TACTed trees")
# )

